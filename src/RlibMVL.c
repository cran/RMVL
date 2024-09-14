#include <stdio.h>
#include <math.h>
#ifndef __WIN32__
#include <sys/mman.h>
#else
#include <windows.h>
#include <winbase.h>
#include <io.h>
#endif
#include <errno.h>
#include "libMVL.h"

/* This needs to be before including R headers because they redefine "error" as a macro */
static inline int mvl_get_error(LIBMVL_CONTEXT *ctx)
{
return ctx->error;
}

#include <R.h>
#include <Rinternals.h>
//#include <Rext/Print.h>

static inline SEXP mkCharU(unsigned char *s)
{
return(mkChar((char *)s));
}

static inline SEXP mkCharLenU(const unsigned char *s, int len)
{
return(mkCharLen((char *)s, len));
}

typedef struct {
	FILE *f;
	char *data;
	LIBMVL_OFFSET64 length;
	LIBMVL_CONTEXT *ctx;
#ifdef __WIN32__
	HANDLE f_handle;
	HANDLE f_map_handle;
#endif
	int modified;
	} MMAPED_LIBRARY;
	
MMAPED_LIBRARY *libraries=NULL;
size_t libraries_size=0;
int libraries_free=0;

SEXP mmap_library(SEXP filename, SEXP mode0)
{
int i, mode;
int idx;
MMAPED_LIBRARY *p;
const char *fn;
SEXP ans;

if(length(mode0)!=1) {
	error("mmap_library argument mode has to be length 1 integer");
	return(R_NilValue);
	}
	
mode=INTEGER(mode0)[0];
fn=CHAR(asChar(filename));

idx=-1;
for(i=0;i<libraries_free;i++) {
	if(libraries[i].ctx==NULL) {
		idx=i;
		break;
		}
	}

if(idx<0 && (libraries_free>=libraries_size)) {
	libraries_size=2*libraries_size+10;
	p=calloc(libraries_size, sizeof(*libraries));
	if(p==NULL) {
		error("Opening MVL library \"%s\": out of memory", fn);
		return(R_NilValue);
		}
	if(libraries_free>0)memcpy(p, libraries, libraries_free*sizeof(*libraries));
	free(libraries);
	libraries=p;
	}
if(idx<0) {
	idx=libraries_free;
	libraries_free++;
	}

//Rprintf("Accessing MVL library from %s\n", fn);

p=&libraries[idx];
memset(p, 0, sizeof(*p));

#ifdef __WIN32__
switch(mode) {
	case 0:
		p->f=fopen(fn, "rb");
		break;
	case 1:
		p->f=fopen(fn, "rb+");
		break;
	case 2: 
		p->f=fopen(fn, "wb");
		break;
	case 3:
		p->f=fopen(fn, "wb+");
		break;
	default:
		error("Unknown mode %d", mode);
		return(R_NilValue);
		
	}
#else 
switch(mode) {
	case 0:
		p->f=fopen(fn, "r");
		break;
	case 1:
		p->f=fopen(fn, "r+");
		break;
	case 2: 
		p->f=fopen(fn, "w");
		break;
	case 3:
		p->f=fopen(fn, "w+");
		break;
	default:
		error("Unknown mode %d", mode);
		return(R_NilValue);
		
	}
#endif
if(p->f==NULL) {
	error("Opening MVL library \"%s\": %s", fn, strerror(errno));
	return(R_NilValue);
	}
	
fseek(p->f, 0, SEEK_END);
p->length=ftell(p->f);
fseek(p->f, 0, SEEK_SET);

p->ctx=mvl_create_context();
p->ctx->f=p->f;

if(p->length>0) {
#ifdef __WIN32__
	p->f_handle=(HANDLE)_get_osfhandle(fileno(p->f));
	if(p->f_handle==NULL) {
		error("Cannot obtain Win32 file handle");
		fclose(p->f);
		p->f=NULL;
		return(R_NilValue);
		}
	
	p->f_map_handle=CreateFileMappingA(p->f_handle, NULL, PAGE_READONLY, 0, 0, NULL);
	if(p->f_map_handle==NULL) {
		error("Cannot obtain Win32 file mapping object");
		fclose(p->f);
		p->f=NULL;
		return(R_NilValue);
		}
		
	p->data=MapViewOfFile(p->f_map_handle, FILE_MAP_READ, 0, 0, p->length);
	if(p->data==NULL) {
		error("Cannot create Win32 File mapping view");
		fclose(p->f);
		p->f=NULL;
		return(R_NilValue);
		}

#else
	p->data=mmap(NULL, p->length, PROT_READ, MAP_SHARED, fileno(p->f), 0);
	if(p->data==NULL) {
		error("Memory mapping MVL library: %s", strerror(errno));
		fclose(p->f);
		p->f=NULL;
		return(R_NilValue);
		}
#endif
	mvl_load_image(p->ctx, p->data, p->length);
	fseek(p->f, 0, SEEK_END);
	if(mode==0) {
		/* Read-only mapping; no need to use up a file descriptor */
		fclose(p->f);
		p->f=NULL;
		p->ctx->f=NULL;
		}
	} else {
	mvl_write_preamble(p->ctx);
	p->modified=1;
	}
	
ans=PROTECT(allocVector(INTSXP, 1));
INTEGER(ans)[0]=idx;
UNPROTECT(1);
return(ans);
}

SEXP remap_library(SEXP idx0, SEXP mode0)
{
int mode;
int idx;
MMAPED_LIBRARY *p;
LIBMVL_OFFSET64 new_length;
size_t cur;

if(length(idx0)!=1) {
	error("close_library requires a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0)return(R_NilValue);
if(idx>=libraries_free)return(R_NilValue);

p=&(libraries[idx]);


if(length(mode0)!=1) {
	error("mmap_library argument mode has to be length 1 integer");
	return(R_NilValue);
	}
	
mode=INTEGER(mode0)[0];

if(p->f==NULL) {
	error("Cannot remap read-only library");
	return(R_NilValue);
	}
	
if(mode==0) {
	/* Read-only mapping; no need to use up a file descriptor */
	if(p->modified) {
		mvl_close(p->ctx);
		if(mvl_get_error(p->ctx)!=0)
			error("Error %d encountered when closing MVL file: %s", mvl_get_error(p->ctx), mvl_strerror(p->ctx));
		}			
	}
	
fflush(p->f);
	
cur=ftell(p->f);
fseek(p->f, 0, SEEK_END);
new_length=ftell(p->f);
fseek(p->f, cur, SEEK_SET);

if(new_length>0) {
#ifdef __WIN32__
	p->f_handle=(HANDLE)_get_osfhandle(fileno(p->f));
	if(p->f_handle==NULL) {
		error("Cannot obtain Win32 file handle");
		fclose(p->f);
		p->f=NULL;
		return(R_NilValue);
		}
	
	UnmapViewOfFile(p->data);
	CloseHandle(p->f_map_handle);
	p->length=new_length;
	
	p->f_map_handle=CreateFileMappingA(p->f_handle, NULL, PAGE_READONLY, 0, 0, NULL);
	if(p->f_map_handle==NULL) {
		error("Cannot obtain Win32 file mapping object");
		fclose(p->f);
		p->f=NULL;
		return(R_NilValue);
		}
		
	p->data=MapViewOfFile(p->f_map_handle, FILE_MAP_READ, 0, 0, p->length);
	if(p->data==NULL) {
		error("Cannot create Win32 File mapping view");
		fclose(p->f);
		p->f=NULL;
		return(R_NilValue);
		}

#else
	if(p->data!=NULL) {
		if(munmap(p->data, p->length)!=0) {
			error("Unmapping data: %s", strerror(errno));
			}
		}

	p->length=new_length;

	p->data=mmap(NULL, p->length, PROT_READ, MAP_SHARED, fileno(p->f), 0);
	if(p->data==NULL) {
		error("Memory mapping MVL library: %s", strerror(errno));
		return(R_NilValue);
		}
#endif
	if(mode==0) {
		/* Read-only mapping; no need to use up a file descriptor */
		fclose(p->f);
		p->f=NULL;
		p->ctx->f=NULL;
		}
	}
	
return(R_NilValue);
}

SEXP close_library(SEXP idx0)
{
int idx; 
MMAPED_LIBRARY *p;
if(length(idx0)!=1) {
	error("close_library requires a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0)return(R_NilValue);
if(idx>=libraries_free)return(R_NilValue);

p=&(libraries[idx]);

if(p->ctx==NULL)return(R_NilValue);

if(p->data!=NULL) {
#ifndef __WIN32__
	if(munmap(p->data, p->length)!=0) {
		error("Unmapping data: %s", strerror(errno));
		}
#else
	UnmapViewOfFile(p->data);
	CloseHandle(p->f_map_handle);
#endif
	p->data=NULL;
	}
	
if(p->modified) {
	mvl_close(p->ctx);
	if(mvl_get_error(p->ctx)!=0)
		error("Error %d encountered when closing MVL file: %s", mvl_get_error(p->ctx), mvl_strerror(p->ctx));
	}
	
mvl_free_context(p->ctx);
p->ctx=NULL;
if(p->f!=NULL)fclose(p->f);
p->f=NULL;

return(R_NilValue);
}


SEXP VECTOR_ELT_STR(SEXP list, const char *s)
{
SEXP elt=R_NilValue;
SEXP names=getAttrib(list, R_NamesSymbol);
SEXP sch;

if(xlength(names)<xlength(list))return(R_NilValue);

for (long i=0; i<xlength(list); i++) {
	sch=STRING_ELT(names, i);
	if(sch==NA_STRING)continue;
	if(!strcmp(CHAR(sch), s)) {
		elt = VECTOR_ELT(list, i);
		break;
		}
	}
return elt;
}


void decode_mvl_object(SEXP obj, int *idx, LIBMVL_OFFSET64 *ofs)
{
SEXP sidx, sofs;
double doffset;
LIBMVL_OFFSET64 *offset=(LIBMVL_OFFSET64 *)&doffset;

sidx=PROTECT(VECTOR_ELT_STR(obj, "handle"));
sofs=VECTOR_ELT_STR(obj, "offset");

*idx=-1;
*ofs=0;

if(sidx!=R_NilValue && length(sidx)==1)*idx=INTEGER(sidx)[0];

if((*idx)>=0 && sofs!=R_NilValue && length(sofs)==1) {
	doffset=REAL(sofs)[0];
	*ofs=*offset;
	}
UNPROTECT(1);
}

LIBMVL_VECTOR * get_mvl_vector(int idx, LIBMVL_OFFSET64 offset)
{
int err;
if(idx<0 || idx>=libraries_free || offset==0)return(NULL);

if(libraries[idx].ctx==NULL)return(NULL);

if(libraries[idx].data==NULL)return(NULL);

if((err=mvl_validate_vector(offset, libraries[idx].data, libraries[idx].length))<0) {
	error("Invalid vector: error %d", err);
	return(NULL);
	}

return((LIBMVL_VECTOR *)(&libraries[idx].data[offset]));
}

LIBMVL_NAMED_LIST * get_mvl_named_list(int idx, LIBMVL_OFFSET64 offset)
{
if(idx<0 || idx>=libraries_free || offset==0)return(NULL);

if(libraries[idx].ctx==NULL)return(NULL);

if(libraries[idx].data==NULL)return(NULL);

return(mvl_read_named_list(libraries[idx].ctx, libraries[idx].data, libraries[idx].length, offset));
}

int get_indices(SEXP indices, LIBMVL_VECTOR *vector, LIBMVL_OFFSET64 *N0, LIBMVL_OFFSET64 **v_idx0)
{
LIBMVL_OFFSET64 *v_idx;
LIBMVL_OFFSET64 N, count, data_offset;
int data_idx;
LIBMVL_VECTOR *vec;

double *pd;
int *pi;

*N0=0;
*v_idx0=NULL;

switch(TYPEOF(indices)) {
	case VECSXP:
		decode_mvl_object(indices, &data_idx, &data_offset);
		vec=get_mvl_vector(data_idx, data_offset);
		
		if(vec==NULL) {
			error("Invalid MVL object or R vector passed as indices");
			return(-1);
			}
		N=mvl_vector_length(vec);
		v_idx=calloc(N, sizeof(*v_idx));
		if(v_idx==NULL) {
			error("Not enough memory");
			return(-2);
			}
		
		switch(mvl_vector_type(vec)) {
			case LIBMVL_VECTOR_OFFSET64:
				for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=mvl_vector_data_int64(vec)[m]-1;
//				memcpy(v_idx, mvl_vector_data_offset(vec), N*sizeof(*v_idx));
				break;
			case LIBMVL_VECTOR_INT32:
				for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=mvl_vector_data_int32(vec)[m]-1;
				break;
			case LIBMVL_VECTOR_INT64:
				for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=mvl_vector_data_int64(vec)[m]-1;
				break;
			case LIBMVL_VECTOR_DOUBLE:
				for(LIBMVL_OFFSET64 m=0;m<N;m++) {
					double a=mvl_vector_data_double(vec)[m];
					v_idx[m]=isfinite(a) ? (LIBMVL_OFFSET64) (a-1.0) : ~(LIBMVL_OFFSET64)0;
					}
				break;
			case LIBMVL_VECTOR_FLOAT:
				for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=mvl_vector_data_float(vec)[m]-1;
				break;
			default:
				error("Cannot interpret MVL object as indices");
				return(-3);
				break;
			}
		break;
	case REALSXP:
		N=xlength(indices);
		v_idx=calloc(N, sizeof(*v_idx));
		if(v_idx==NULL) {
			error("Not enough memory");
			return(-4);
			}
		pd=REAL(indices);
		for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=pd[m]-1;
		break;
	case INTSXP:
		N=xlength(indices);
		v_idx=calloc(N, sizeof(*v_idx));
		if(v_idx==NULL) {
			error("Not enough memory");
			return(-5);
			}
		pi=INTEGER(indices);
		for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=pi[m]-1;
		break;
	case LGLSXP:
		N=xlength(indices);
		pi=LOGICAL(indices);
		count=0;
		for(LIBMVL_OFFSET64 m=0;m<N;m++)if(pi[m])count++;
		
		v_idx=calloc(count, sizeof(*v_idx));
		if(v_idx==NULL) {
			error("Not enough memory");
			return(-5);
			}
		for(LIBMVL_OFFSET64 m=0, k=0;m<N;m++)
			if(pi[m] && (pi[m]!=NA_LOGICAL)) {
				v_idx[k]=m;
				k++;
				}
		break;
	case NILSXP:
		if(vector==NULL) {
			error("Cannot infer vector length");
			return(-6);
			}
		N=mvl_vector_length(vector);
		if(mvl_vector_type(vector)==LIBMVL_PACKED_LIST64)N--;
		v_idx=calloc(N, sizeof(*v_idx));
		if(v_idx==NULL) {
			error("Not enough memory");
			return(-7);
			}
		for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=m;
		break;
	default:
		error("Cannot interpret R object as index");
		return(-8);		
	}
*N0=N;
*v_idx0=v_idx;
return(0);
}

SEXP get_status(void)
{
int max_N=10;
int idx=0, idx2, n_open, nunprot=0;
int j;
SEXP names, ans, data;
	
names=PROTECT(allocVector(STRSXP, max_N));
ans=PROTECT(allocVector(VECSXP, max_N));

#define add_status(name, val)   {\
	if(idx<max_N) { \
		SET_VECTOR_ELT(ans, idx, PROTECT(val)); \
		SET_STRING_ELT(names, idx, mkChar(name)); \
		idx++; \
		} else { \
		Rprintf("*** MVL_INTERNAL_ERROR: too many status fields\n"); \
		} \
	}

add_status("size_t_bytes", ScalarInteger(sizeof(size_t)));
add_status("off_t_bytes", ScalarInteger(sizeof(off_t)));
add_status("long_bytes", ScalarInteger(sizeof(long)));
add_status("offset64_bytes", ScalarInteger(sizeof(LIBMVL_OFFSET64)));
add_status("vector_header_bytes", ScalarInteger(sizeof(LIBMVL_VECTOR_HEADER)));

// This confuses Rchk, but doing manually is worse because Rchk does not properly handle macro expansion

n_open=0;
for(j=0;j<libraries_free;j++)
	if(libraries[j].ctx!=NULL)n_open++;
	
add_status("open_libraries", ScalarInteger(n_open));
	
data=PROTECT(allocVector(INTSXP, n_open));
idx2=0;
for(j=0;j<libraries_free;j++)
	if(libraries[j].ctx!=NULL) {
		INTEGER(data)[idx2]=j;
		idx2++;
		}
		
add_status("library_handles", data);
nunprot++;

data=PROTECT(allocVector(INTSXP, n_open));
idx2=0;
for(j=0;j<libraries_free;j++)
	if(libraries[j].ctx!=NULL) {
		INTEGER(data)[idx2]=libraries[j].ctx->flags;
		idx2++;
		}
		
add_status("library_flags", data);
nunprot++;

data=PROTECT(allocVector(LGLSXP, n_open));
idx2=0;
for(j=0;j<libraries_free;j++)
	if(libraries[j].ctx!=NULL) {
		LOGICAL(data)[idx2]=libraries[j].modified;
		idx2++;
		}
		
add_status("library_modified", data);
nunprot++;

data=PROTECT(allocVector(REALSXP, n_open));
idx2=0;
for(j=0;j<libraries_free;j++)
	if(libraries[j].ctx!=NULL) {
		REAL(data)[idx2]=libraries[j].length;
		idx2++;
		}
		
add_status("library_length", data);
nunprot++;

#undef add_status
if(idx!=max_N) {
	Rprintf("*** RMVL INTERNAL ERROR: idx=%d vs max_N=%d in %s:%d\n",
		idx, max_N, __FILE__, __LINE__);
	}

setAttrib(ans, R_NamesSymbol, names);
UNPROTECT(idx+nunprot+2);
return(ans);
}

SEXP find_directory_entries(SEXP idx0, SEXP tag)
{
int idx;
SEXP ans, class, sch;
const char *tag0;
long i;
LIBMVL_OFFSET64 offset;
double *doffset=(double *)&offset;
if(length(idx0)!=1) {
	error("find_directory_entry first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL){
	error("invalid MVL handle");
	return(R_NilValue);
	}
ans=PROTECT(allocVector(REALSXP, xlength(tag)));
for(i=0;i<xlength(tag);i++) {
	sch=STRING_ELT(tag, i);
	if(sch==NA_STRING) {
		offset=0;
		REAL(ans)[i]=*doffset;
		} else {
		tag0=CHAR(sch);
		offset=mvl_find_directory_entry(libraries[idx].ctx, tag0);
		REAL(ans)[i]=*doffset;
		}
	}

class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(2);
return(ans);
}

SEXP get_directory(SEXP idx0)
{
int idx;
SEXP ans, class, names;
long i;
LIBMVL_OFFSET64 offset;
double *doffset=(double *)&offset;
LIBMVL_NAMED_LIST *dir;

double *dp;

if(length(idx0)!=1) {
	error("find_directory_entry first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL){
	error("invalid MVL handle");
	return(R_NilValue);
	}
dir=libraries[idx].ctx->directory;
ans=PROTECT(allocVector(REALSXP, dir->free));
names=PROTECT(allocVector(STRSXP, dir->free));
dp=REAL(ans);
for(i=0;i<dir->free;i++) {
	SET_STRING_ELT(names, i, mkCharLenU(dir->tag[i], dir->tag_length[i]));
	offset=dir->offset[i];
	dp[i]=*doffset;
	}
setAttrib(ans, R_NamesSymbol, names);
class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(3);
return(ans);
}

SEXP read_lengths(SEXP idx0, SEXP offsets)
{
int idx;
SEXP ans;
long i;
double doffset;
LIBMVL_OFFSET64 *offset0=(LIBMVL_OFFSET64 *)&doffset;
LIBMVL_OFFSET64 offset;
LIBMVL_VECTOR *vec;
double *d_ans, *d_offsets;

if(length(idx0)!=1) {
	error("find_directory_entry first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
ans=PROTECT(allocVector(REALSXP, xlength(offsets)));
d_ans=REAL(ans);
d_offsets=REAL(offsets);
for(i=0;i<xlength(offsets);i++) {
	doffset=d_offsets[i];
	offset=*offset0;
	if(mvl_validate_vector(offset, libraries[idx].data, libraries[idx].length)) {
		d_ans[i]=NA_REAL;
		continue;
		}
	vec=(LIBMVL_VECTOR *)(&libraries[idx].data[offset]);
	switch(mvl_vector_type(vec)) {
		case LIBMVL_PACKED_LIST64:
			d_ans[i]=mvl_vector_length(vec)-1;
			break;
		default:
			d_ans[i]=mvl_vector_length(vec);
			break;
		}
	}

UNPROTECT(1);
return(ans);
}

SEXP compute_vector_stats(SEXP idx0, SEXP offsets)
{
int idx;
SEXP ans;
long i, j;
double doffset;
LIBMVL_OFFSET64 *offset0=(LIBMVL_OFFSET64 *)&doffset;
LIBMVL_OFFSET64 offset;
LIBMVL_VECTOR *vec;
double *d_ans, *d_offsets;
LIBMVL_VEC_STATS stats;
int nfields=sizeof(stats)/sizeof(double);

if(length(idx0)!=1) {
	error("find_directory_entry first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
ans=PROTECT(allocVector(REALSXP, xlength(offsets)*nfields));
d_ans=REAL(ans);
d_offsets=REAL(offsets);
for(i=0;i<xlength(offsets);i++) {
	doffset=d_offsets[i];
	offset=*offset0;
	if(mvl_validate_vector(offset, libraries[idx].data, libraries[idx].length)) {		
		for(j=0;j<nfields;j++)
			d_ans[i*nfields+j]=NA_REAL;
		continue;
		}
	vec=(LIBMVL_VECTOR *)(&libraries[idx].data[offset]);
	mvl_compute_vec_stats(vec, &stats);
	memcpy(&(d_ans[i*nfields]), &stats, sizeof(stats));
	}

UNPROTECT(1);
return(ans);
}


SEXP read_types(SEXP idx0, SEXP offsets)
{
int idx;
SEXP ans;
long i;
double doffset;
LIBMVL_OFFSET64 *offset0=(LIBMVL_OFFSET64 *)&doffset;
LIBMVL_OFFSET64 offset;
LIBMVL_VECTOR *vec;

int *p_ans;
double *p_offsets;

if(length(idx0)!=1) {
	error("find_directory_entry first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
ans=PROTECT(allocVector(INTSXP, xlength(offsets)));
p_ans=INTEGER(ans);
p_offsets=REAL(offsets);
for(i=0;i<xlength(offsets);i++) {
	doffset=p_offsets[i];
	offset=*offset0;
	if(mvl_validate_vector(offset, libraries[idx].data, libraries[idx].length)) {		
		p_ans[i]=NA_INTEGER;
		continue;
		}
	vec=(LIBMVL_VECTOR *)(&libraries[idx].data[offset]);
	p_ans[i]=mvl_vector_type(vec);
	}

UNPROTECT(1);
return(ans);
}

/* Prefer raw vectors for types with no exact R equivalent */
SEXP read_vectors_raw(SEXP idx0, SEXP offsets)
{
int idx;
SEXP ans, v, class;
long i, j;
long field_size;
double doffset;
LIBMVL_OFFSET64 *offset0=(LIBMVL_OFFSET64 *)&doffset;
LIBMVL_OFFSET64 offset;
double *doffset2=(double *)&offset;
LIBMVL_VECTOR *vec;

double *pd;
int *pi;
unsigned char *pc;

if(length(idx0)!=1) {
	error("find_directory_entry first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
ans=PROTECT(allocVector(VECSXP, xlength(offsets)));
for(i=0;i<xlength(offsets);i++) {
	doffset=REAL(offsets)[i];
	offset=*offset0;
	if(offset==0 || offset>libraries[idx].length-sizeof(LIBMVL_VECTOR_HEADER)) {
		SET_VECTOR_ELT(ans, i, R_NilValue);
		continue;
		}
	vec=(LIBMVL_VECTOR *)(&libraries[idx].data[offset]);
	switch(mvl_vector_type(vec)) {
		case LIBMVL_VECTOR_UINT8:
		case LIBMVL_VECTOR_CSTRING:
				field_size=1;
				break;
		case LIBMVL_VECTOR_FLOAT:
		case LIBMVL_VECTOR_INT32:
				field_size=4;
				break;
		case LIBMVL_VECTOR_DOUBLE:
		case LIBMVL_VECTOR_INT64:
		case LIBMVL_VECTOR_OFFSET64:
				field_size=8;
				break;
		default:
			field_size=1;
		}
	switch(mvl_vector_type(vec)) {
		case LIBMVL_VECTOR_UINT8:
		case LIBMVL_VECTOR_INT64:
		case LIBMVL_VECTOR_FLOAT:
			v=PROTECT(allocVector(RAWSXP, mvl_vector_length(vec)*field_size));
			pc=RAW(v);
			for(j=0;j<mvl_vector_length(vec)*field_size;j++)
				pc[j]=mvl_vector_data_uint8(vec)[j];
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_CSTRING:
			v=PROTECT(allocVector(STRSXP, 1));
			/* TODO: check that vector length is within R limits */
			if(mvl_string_is_na((char *)mvl_vector_data_uint8(vec), mvl_vector_length(vec)))
				SET_STRING_ELT(v, 0, NA_STRING);
				else
				SET_STRING_ELT(v, 0, mkCharLenU(mvl_vector_data_uint8(vec), mvl_vector_length(vec)));
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			//SET_VECTOR_ELT(ans, i, mkCharLenU(mvl_vector_data_uint8(v), mvl_vector_length(vec)));
			break;
		case LIBMVL_VECTOR_INT32:
			v=PROTECT(allocVector(INTSXP, mvl_vector_length(vec)));
			pi=INTEGER(v);
			for(j=0;j<mvl_vector_length(vec);j++)
				pi[j]=mvl_vector_data_int32(vec)[j];
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_DOUBLE:
			v=PROTECT(allocVector(REALSXP, mvl_vector_length(vec)));
			pd=REAL(v);
			for(j=0;j<mvl_vector_length(vec);j++)
				pd[j]=mvl_vector_data_double(vec)[j];
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_OFFSET64:
			v=PROTECT(allocVector(REALSXP, mvl_vector_length(vec)));
			pd=REAL(v);
			for(j=0;j<mvl_vector_length(vec);j++) {
				offset=mvl_vector_data_offset(vec)[j];
				pd[j]=*doffset2;
				}
			class=PROTECT(allocVector(STRSXP, 1));
			SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
			classgets(v, class);
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(2);
			break;
		case LIBMVL_PACKED_LIST64:
			v=PROTECT(allocVector(STRSXP, mvl_vector_length(vec)-1));
			/* TODO: check that vector length is within R limits */
			for(j=0;j<mvl_vector_length(vec)-1;j++) {
				if(mvl_packed_list_is_na(vec, libraries[idx].data, j))
					SET_STRING_ELT(v, j, NA_STRING);
					else
					SET_STRING_ELT(v, j, mkCharLenU(mvl_packed_list_get_entry(vec, libraries[idx].data, j), mvl_packed_list_get_entry_bytelength(vec, j)));
				}
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		default:
			warning("Unknown vector type");
			SET_VECTOR_ELT(ans, i, R_NilValue);
			break;
		}
		
	}

UNPROTECT(1);
return(ans);
}

SEXP read_vectors_idx_raw2(SEXP idx0, SEXP offsets, SEXP indicies)
{
int idx;
SEXP ans, v, class;
long i, j, field_size;
double doffset;
LIBMVL_OFFSET64 *offset0=(LIBMVL_OFFSET64 *)&doffset;
LIBMVL_OFFSET64 offset;
double *doffset2=(double *)&offset;
LIBMVL_VECTOR *vec;
LIBMVL_OFFSET64 *v_idx, N, m, k;

double *pd;
int *pi;
unsigned char *pc;

if(length(idx0)!=1) {
	error("read_vectors_idx_raw2 first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
	
if(get_indices(indicies, NULL, &N, &v_idx)!=0) {
	return(R_NilValue);
	}
	
ans=PROTECT(allocVector(VECSXP, xlength(offsets)));
for(i=0;i<xlength(offsets);i++) {
	doffset=REAL(offsets)[i];
	offset=*offset0;
	if(offset==0 || offset>libraries[idx].length-sizeof(LIBMVL_VECTOR_HEADER)) {
		SET_VECTOR_ELT(ans, i, R_NilValue);
		continue;
		}
	vec=(LIBMVL_VECTOR *)(&libraries[idx].data[offset]);

	m=mvl_vector_length(vec);
	for(j=0;j<N;j++) {
		k=v_idx[j];
		if((k<0) || (k>m)) {
			UNPROTECT(1);
			error("Index is out of range");
			free(v_idx);
			return(R_NilValue);
			}
		}
	
	field_size=mvl_element_size(mvl_vector_type(vec));
	
	/* Unknown type */
	if(field_size<1) {
		SET_VECTOR_ELT(ans, i, R_NilValue);
		continue;
		}

	switch(mvl_vector_type(vec)) {
		case LIBMVL_VECTOR_UINT8:
			v=PROTECT(allocVector(RAWSXP, N));
			pc=RAW(v);
			for(j=0;j<N;j++)
				pc[j]=mvl_vector_data_uint8(vec)[v_idx[j]];
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			//SET_VECTOR_ELT(ans, i, mkCharLenU(mvl_vector_data_uint8(v), mvl_vector_length(vec)));
			break;
		case LIBMVL_VECTOR_CSTRING:
			error("String subset not supported");
			free(v_idx);
			return(R_NilValue);
			v=PROTECT(allocVector(STRSXP, 1));
			/* TODO: check that vector length is within R limits */
			SET_STRING_ELT(v, 0, mkCharLenU(mvl_vector_data_uint8(vec), mvl_vector_length(vec)));
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			//SET_VECTOR_ELT(ans, i, mkCharLenU(mvl_vector_data_uint8(v), mvl_vector_length(vec)));
			break;
		case LIBMVL_VECTOR_INT32:
			v=PROTECT(allocVector(INTSXP, N));
			pi=INTEGER(v);
			for(j=0;j<N;j++)
				pi[j]=mvl_vector_data_int32(vec)[v_idx[j]];
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_INT64:
			v=PROTECT(allocVector(RAWSXP, N*field_size));
			pc=RAW(v);
			for(j=0;j<N*field_size;j+=field_size)
				memcpy(&(pc[j]), &(mvl_vector_data_int64(vec)[v_idx[j]]), field_size);
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);			
			break;
		case LIBMVL_VECTOR_FLOAT:
			v=PROTECT(allocVector(RAWSXP, N*field_size));
			pc=RAW(v);
			for(j=0;j<N*field_size;j+=field_size)
				memcpy(&(pc[j]), &(mvl_vector_data_float(vec)[v_idx[j]]), field_size);
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);			
			break;
		case LIBMVL_VECTOR_DOUBLE:
			v=PROTECT(allocVector(REALSXP, N));
			pd=REAL(v);
			for(j=0;j<N;j++)
				pd[j]=mvl_vector_data_double(vec)[v_idx[j]];
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_OFFSET64:
			v=PROTECT(allocVector(REALSXP, N));
			pd=REAL(v);
			for(j=0;j<N;j++) {
				offset=mvl_vector_data_offset(vec)[v_idx[j]];
				pd[j]=*doffset2;
				}
			class=PROTECT(allocVector(STRSXP, 1));
			SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
			classgets(v, class);
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(2);
			break;
		case LIBMVL_PACKED_LIST64:
			v=PROTECT(allocVector(STRSXP, N));
			/* TODO: check that vector length is within R limits */
			for(j=0;j<N;j++) {
				if(mvl_packed_list_is_na(vec, libraries[idx].data, v_idx[j]))
					SET_STRING_ELT(v, j, NA_STRING);
					else
					SET_STRING_ELT(v, j, mkCharLenU(mvl_packed_list_get_entry(vec, libraries[idx].data, v_idx[j]), mvl_packed_list_get_entry_bytelength(vec, v_idx[j])));
				}
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		default:
			warning("Unknown vector type");
			SET_VECTOR_ELT(ans, i, R_NilValue);
			break;
		}
	}

UNPROTECT(1);
free(v_idx);
return(ans);
}


SEXP get_vector_data_ptr(SEXP idx0, SEXP offsets)
{
int idx;
SEXP ans;
long i;
double doffset;
LIBMVL_OFFSET64 *offset0=(LIBMVL_OFFSET64 *)&doffset;
LIBMVL_OFFSET64 offset;
double *doffset2=(double *)&offset;
LIBMVL_VECTOR *vec;

double *p_offsets, *p_ans;

if(length(idx0)!=1) {
	error("read_vectors_idx_raw2 first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
ans=PROTECT(allocVector(REALSXP, xlength(offsets)));
p_ans=REAL(ans);
p_offsets=REAL(offsets);
for(i=0;i<xlength(offsets);i++) {
	doffset=p_offsets[i];
	offset=*offset0;
	if(mvl_validate_vector(offset, libraries[idx].data, libraries[idx].length)) {
		UNPROTECT(1);
		error("Invalid vector offset provided");
		return(R_NilValue);
		}
	vec=(LIBMVL_VECTOR *)(&libraries[idx].data[offset]);
	
	offset=(LIBMVL_OFFSET64)&(mvl_vector_data(vec));
	p_ans[i]=*doffset2;	
	}
UNPROTECT(1);
return(ans);
}

SEXP read_vectors(SEXP idx0, SEXP offsets)
{
int idx;
SEXP ans, v, class;
long i, j;
double doffset;
LIBMVL_OFFSET64 *offset0=(LIBMVL_OFFSET64 *)&doffset;
LIBMVL_OFFSET64 offset;
double *doffset2=(double *)&offset;
LIBMVL_VECTOR *vec;

double *pd;
int *pi;
unsigned char *pc;

if(length(idx0)!=1) {
	error("read_vectors first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
ans=PROTECT(allocVector(VECSXP, xlength(offsets)));
for(i=0;i<xlength(offsets);i++) {
	doffset=REAL(offsets)[i];
	offset=*offset0;
	if(offset==0 || offset>libraries[idx].length-sizeof(LIBMVL_VECTOR_HEADER)) {
		SET_VECTOR_ELT(ans, i, R_NilValue);
		continue;
		}
	vec=(LIBMVL_VECTOR *)(&libraries[idx].data[offset]);
	switch(mvl_vector_type(vec)) {
		case LIBMVL_VECTOR_UINT8:
			v=PROTECT(allocVector(RAWSXP, mvl_vector_length(vec)));
			pc=RAW(v);
			for(j=0;j<mvl_vector_length(vec);j++)
				pc[j]=mvl_vector_data_uint8(vec)[j];
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			//SET_VECTOR_ELT(ans, i, mkCharLenU(mvl_vector_data_uint8(v), mvl_vector_length(vec)));
			break;
		case LIBMVL_VECTOR_CSTRING:
			v=PROTECT(allocVector(STRSXP, 1));
			/* TODO: check that vector length is within R limits */
			SET_STRING_ELT(v, 0, mkCharLenU(mvl_vector_data_uint8(vec), mvl_vector_length(vec)));
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			//SET_VECTOR_ELT(ans, i, mkCharLenU(mvl_vector_data_uint8(v), mvl_vector_length(vec)));
			break;
		case LIBMVL_VECTOR_INT32:
			v=PROTECT(allocVector(INTSXP, mvl_vector_length(vec)));
			pi=INTEGER(v);
			for(j=0;j<mvl_vector_length(vec);j++)
				pi[j]=mvl_vector_data_int32(vec)[j];
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_INT64:
			warning("Converted 64-bit integers to doubles");
			v=PROTECT(allocVector(REALSXP, mvl_vector_length(vec)));
			pd=REAL(v);
			for(j=0;j<mvl_vector_length(vec);j++)
				pd[j]=mvl_vector_data_int64(vec)[j];
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_FLOAT:
			//warning("Converted 32-bit floats to doubles");
			v=PROTECT(allocVector(REALSXP, mvl_vector_length(vec)));
			pd=REAL(v);
			for(j=0;j<mvl_vector_length(vec);j++)
				pd[j]=mvl_vector_data_float(vec)[j];
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_DOUBLE:
			v=PROTECT(allocVector(REALSXP, mvl_vector_length(vec)));
			pd=REAL(v);
			for(j=0;j<mvl_vector_length(vec);j++)
				pd[j]=mvl_vector_data_double(vec)[j];
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_OFFSET64:
			v=PROTECT(allocVector(REALSXP, mvl_vector_length(vec)));
			pd=REAL(v);
			for(j=0;j<mvl_vector_length(vec);j++) {
				offset=mvl_vector_data_offset(vec)[j];
				pd[j]=*doffset2;
				}
			class=PROTECT(allocVector(STRSXP, 1));
			SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
			classgets(v, class);
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(2);
			break;
		case LIBMVL_PACKED_LIST64:
			v=PROTECT(allocVector(STRSXP, mvl_vector_length(vec)-1));
			/* TODO: check that vector length is within R limits */
			for(j=0;j<mvl_vector_length(vec)-1;j++) {
				if(mvl_packed_list_is_na(vec, libraries[idx].data, j))
					SET_STRING_ELT(v, j, NA_STRING);
					else
					SET_STRING_ELT(v, j, mkCharLenU(mvl_packed_list_get_entry(vec, libraries[idx].data, j), mvl_packed_list_get_entry_bytelength(vec, j)));
				}
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		default:
			warning("Unknown vector type");
			SET_VECTOR_ELT(ans, i, R_NilValue);
			break;
		}
	}

UNPROTECT(1);
return(ans);
}

static inline LIBMVL_OFFSET64 indexR2c(LIBMVL_OFFSET64 idx, LIBMVL_OFFSET64 idx_default)
{
return(idx<1 ? idx_default : idx-1);
}

/* This function accepts vectors in variety of formats */
SEXP read_vectors_idx3(SEXP idx0, SEXP offsets, SEXP indicies)
{
int idx;
SEXP ans, v, class;
long i;
double doffset;
LIBMVL_OFFSET64 *offset0=(LIBMVL_OFFSET64 *)&doffset;
LIBMVL_OFFSET64 offset;
LIBMVL_VECTOR *vec, *vec_idx=NULL;
LIBMVL_OFFSET64 N, N0;

double * restrict pd;
int * restrict pi;
unsigned char * restrict pc;
LIBMVL_OFFSET64 * restrict poffs;
double * restrict pd2;

if(length(idx0)!=1) {
	error("read_vectors_idx3 first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
	
if(TYPEOF(indicies)==VECSXP) {
	int data_idx;
	LIBMVL_OFFSET64 data_offset;
	decode_mvl_object(indicies, &data_idx, &data_offset);
	vec_idx=get_mvl_vector(data_idx, data_offset);
		
	if(vec_idx==NULL) {
		error("Invalid MVL object or R vector passed as indices");
		return(R_NilValue);
		}
	}
	

ans=PROTECT(allocVector(VECSXP, xlength(offsets)));
for(i=0;i<xlength(offsets);i++) {
	doffset=REAL(offsets)[i];
	offset=*offset0;
	if(offset==0 || offset>libraries[idx].length-sizeof(LIBMVL_VECTOR_HEADER)) {
		SET_VECTOR_ELT(ans, i, R_NilValue);
		continue;
		}
	vec=(LIBMVL_VECTOR *)(&libraries[idx].data[offset]);
	N0=mvl_vector_length(vec);
	
		
	switch(TYPEOF(indicies)) { 
		case VECSXP: {
			
			switch(mvl_vector_type(vec_idx)) {
				case LIBMVL_VECTOR_OFFSET64:
				case LIBMVL_VECTOR_INT32:
				case LIBMVL_VECTOR_INT64:
				case LIBMVL_VECTOR_DOUBLE:
				case LIBMVL_VECTOR_FLOAT:
					N=mvl_vector_length(vec_idx);
					break;
				default:
					error("Cannot interpret MVL object as indices");
					UNPROTECT(1);
					return(R_NilValue);
					break;
				}
			
			break;
			}
		case REALSXP:
		case INTSXP: {
			N=xlength(indicies);
			break; 
			} 
		case LGLSXP: {
			N=0;
			int *pi=LOGICAL(indicies);
			for(LIBMVL_OFFSET64 j=0;j<xlength(indicies);j++)
				if(pi[j] && (pi[j]!=NA_LOGICAL))N++;
			break; 
			}
		case NILSXP: { 
			N=N0;
			break; 
			} 
		default:
			N=0;
			error("Cannot handle R index type %d", TYPEOF(indicies));
			UNPROTECT(1);
			return(R_NilValue);
			break; 
		}
		
#define INDEX_LOOP(line, line_default) \
	{ \
	switch(TYPEOF(indicies)) { \
		case VECSXP: { \
			switch(mvl_vector_type(vec_idx)) { \
				case LIBMVL_VECTOR_OFFSET64: \
					for(LIBMVL_OFFSET64 j=0;j<N;j++) { \
						LIBMVL_OFFSET64 j0=indexR2c(mvl_vector_data_int64(vec_idx)[j], N0); \
						if(j0<N0) { \
							line ;\
							} else { \
							line_default ;\
							} \
						} \
					break; \
				case LIBMVL_VECTOR_INT32: \
					for(LIBMVL_OFFSET64 j=0;j<N;j++) { \
						LIBMVL_OFFSET64 j0=indexR2c(mvl_vector_data_int32(vec_idx)[j], N0); \
						if(j0<N0) { \
							line ;\
							} else { \
							line_default ;\
							} \
						} \
					break; \
				case LIBMVL_VECTOR_INT64: \
					for(LIBMVL_OFFSET64 j=0;j<N;j++) { \
						LIBMVL_OFFSET64 j0=indexR2c(mvl_vector_data_int64(vec_idx)[j], N0); \
						if(j0<N0) { \
							line ;\
							} else { \
							line_default ;\
							} \
						} \
					break; \
				case LIBMVL_VECTOR_DOUBLE: \
					for(LIBMVL_OFFSET64 j=0;j<N;j++) { \
						double vpidx=mvl_vector_data_double(vec_idx)[j]; \
						LIBMVL_OFFSET64 j0=isfinite(vpidx) ? indexR2c(vpidx, N0) : N0; \
						if(j0<N0) { \
							line ;\
							} else { \
							line_default ;\
							} \
						} \
					break; \
				case LIBMVL_VECTOR_FLOAT: \
					for(LIBMVL_OFFSET64 j=0;j<N;j++) { \
						LIBMVL_OFFSET64 j0=indexR2c(mvl_vector_data_float(vec_idx)[j], N0); \
						if(j0<N0) { \
							line ;\
							} else { \
							line_default ;\
							} \
						} \
					break; \
				default: \
					break; \
				} \
			break; \
			} \
		case REALSXP: {\
			double * restrict pidx=REAL(indicies); \
			for(LIBMVL_OFFSET64 j=0;j<N;j++) { \
				double vpidx=pidx[j]; \
				LIBMVL_OFFSET64 j0=isfinite(vpidx) ? indexR2c(vpidx, N0) : N0; \
				if(j0<N0) { \
					line ;\
					} else { \
					line_default ;\
					} \
				} \
			break; \
			} \
		case INTSXP: { \
			int * restrict pidx=INTEGER(indicies); \
			for(LIBMVL_OFFSET64 j=0;j<N;j++) { \
				LIBMVL_OFFSET64 j0=indexR2c(pidx[j], N0); \
				if(j0<N0) { \
					line ;\
					} else { \
					line_default ;\
					} \
				} \
			break; \
			} \
		case LGLSXP: { \
			int * restrict pidx=LOGICAL(indicies); \
			for(LIBMVL_OFFSET64 j0=0,j=0;j0<xlength(indicies);j0++) \
				if(pidx[j0] && (pidx[j0]!=NA_LOGICAL)) { \
					line; \
					j++; \
					} else { \
					line_default ;\
					} \
			break; \
			} \
		case NILSXP: { \
			for(LIBMVL_OFFSET64 j=0;j<N;j++) { \
				LIBMVL_OFFSET64 j0=j; \
					{ \
					line ;\
					} \
				} \
			break; \
			} \
		default:\
			break; \
		} \
	}

	
	switch(mvl_vector_type(vec)) {
		case LIBMVL_VECTOR_UINT8:
			v=PROTECT(allocVector(RAWSXP, N));
			pc=RAW(v);
			INDEX_LOOP(pc[j]=mvl_vector_data_uint8(vec)[j0], pc[j]=0);
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			//SET_VECTOR_ELT(ans, i, mkCharLenU(mvl_vector_data_uint8(v), mvl_vector_length(vec)));
			break;
		case LIBMVL_VECTOR_CSTRING:
			error("String subset not supported");
			return(R_NilValue);
			v=PROTECT(allocVector(STRSXP, 1));
			/* TODO: check that vector length is within R limits */
			SET_STRING_ELT(v, 0, mkCharLenU(mvl_vector_data_uint8(vec), mvl_vector_length(vec)));
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			//SET_VECTOR_ELT(ans, i, mkCharLenU(mvl_vector_data_uint8(v), mvl_vector_length(vec)));
			break;
		case LIBMVL_VECTOR_INT32:
			v=PROTECT(allocVector(INTSXP, N));
			pi=INTEGER(v);
			INDEX_LOOP(pi[j]=mvl_vector_data_int32(vec)[j0], pi[j]=NA_INTEGER);
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_INT64:
			warning("Converted 64-bit integers to doubles");
			v=PROTECT(allocVector(REALSXP, N));
			pd=REAL(v);
			INDEX_LOOP(pd[j]=mvl_vector_data_int64(vec)[j0], pd[j]=NA_REAL);
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_FLOAT:
			//warning("Converted 32-bit floats to doubles");
			v=PROTECT(allocVector(REALSXP, N));
			pd=REAL(v);
			INDEX_LOOP(pd[j]=mvl_vector_data_float(vec)[j0], pd[j]=NA_REAL);
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_DOUBLE:
			v=PROTECT(allocVector(REALSXP, N));
			pd=REAL(v);
			pd2=mvl_vector_data_double(vec);
			INDEX_LOOP(pd[j]=pd2[j0], pd[j]=NA_REAL);
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		case LIBMVL_VECTOR_OFFSET64:
			v=PROTECT(allocVector(REALSXP, N));
			poffs=(LIBMVL_OFFSET64 *)REAL(v);
			INDEX_LOOP(poffs[j]=mvl_vector_data_offset(vec)[j0], poffs[j]=0);
			class=PROTECT(allocVector(STRSXP, 1));
			SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
			classgets(v, class);
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(2);
			break;
		case LIBMVL_PACKED_LIST64:
			v=PROTECT(allocVector(STRSXP, N));
			/* TODO: check that vector length is within R limits */
			INDEX_LOOP( {
				void *data=libraries[idx].data;
				if(mvl_packed_list_is_na(vec, data, j0))
					SET_STRING_ELT(v, j, NA_STRING);
					else
					SET_STRING_ELT(v, j, mkCharLenU(mvl_packed_list_get_entry(vec, data, j0), mvl_packed_list_get_entry_bytelength(vec, j0)));
				}, SET_STRING_ELT(v, j, NA_STRING));
			SET_VECTOR_ELT(ans, i, v);
			UNPROTECT(1);
			break;
		default:
			warning("Unknown vector type");
			SET_VECTOR_ELT(ans, i, R_NilValue);
			break;
		}
	}
#undef INDEX_LOOP
	
UNPROTECT(1);
return(ans);
}


SEXP read_metadata(SEXP idx0, SEXP offsets)
{
int idx;
SEXP ans, class;
long i;
double doffset;
LIBMVL_OFFSET64 *offset0=(LIBMVL_OFFSET64 *)&doffset;
LIBMVL_OFFSET64 offset;
double *doffset2=(double *)&offset;
LIBMVL_VECTOR *vec;
double *dans;
if(length(idx0)!=1) {
	error("read_metadata first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
ans=PROTECT(allocVector(REALSXP, xlength(offsets)));
dans=REAL(ans);
for(i=0;i<xlength(offsets);i++) {
	doffset=REAL(offsets)[i];
	offset=*offset0;
	if(mvl_validate_vector(offset, libraries[idx].data, libraries[idx].length)) {
		Rprintf("offset=%lld data=%p length=%lld\n", offset, libraries[idx].data, libraries[idx].length);
		offset=0;
		dans[i]=NA_REAL;
		continue;
		} 
	vec=(LIBMVL_VECTOR *)(&libraries[idx].data[offset]);
	offset=mvl_vector_metadata_offset(vec);
	dans[i]=*doffset2;
	}

class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(2);
return(ans);
}

SEXP add_directory_entries(SEXP idx0, SEXP tags, SEXP offsets)
{
long i;
int idx;
double doffset;
SEXP sch;
LIBMVL_OFFSET64 *offset=(LIBMVL_OFFSET64 *)&doffset;
if(length(idx0)!=1) {
	error("add_directory_entries first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(xlength(tags)!=xlength(offsets)) {
	error("add_directory_entries requires number of tags to match number of offsets");
	return(R_NilValue);
	}
for(i=0;i<xlength(tags);i++) {
	doffset=REAL(offsets)[i];
	sch=STRING_ELT(tags, i);
	if(sch==NA_STRING)
		warning("Ignoring attempt to add directory entry with NA (missing value) tag");
		else
		mvl_add_directory_entry(libraries[idx].ctx, *offset, CHAR(sch));
	}
return(R_NilValue);
}

SEXP write_vector(SEXP idx0, SEXP type0, SEXP data, SEXP metadata_offset)
{
long i;
int idx, type;
double dmoffset;
LIBMVL_OFFSET64 *moffset=(LIBMVL_OFFSET64 *)&dmoffset;

LIBMVL_OFFSET64 offset;
double *doffset=(double *)&offset;
const char *ch, **strvec2;
LIBMVL_OFFSET64 *strvec;
long long *idata;
long *strlen2;
float *fdata;
unsigned char *bdata;

double *pd;
int *pi;

SEXP ans, class;

if(length(idx0)!=1) {
	error("write_vector first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].f==NULL) {
	error("library not open for writing");
	return(R_NilValue);
	}
if(mvl_get_error(libraries[idx].ctx)!=0) {
	error("library has error status %d: %s", mvl_get_error(libraries[idx].ctx), mvl_strerror(libraries[idx].ctx));
	return(R_NilValue);
	}

if(length(type0)!=1) {
	error("write_vector second argument must be a single integer");
	return(R_NilValue);
	}
type=INTEGER(type0)[0];

libraries[idx].modified=1;

if(length(metadata_offset)<1) {
	*moffset=0;
	} else {
	dmoffset=REAL(metadata_offset)[0];
	}
switch(type) {
	case LIBMVL_VECTOR_UINT8:
		switch(TYPEOF(data)) {
			case RAWSXP:
				offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, xlength(data), RAW(data), *moffset);
				break;
			case LGLSXP:
				bdata=calloc(xlength(data), 1);
				if(bdata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				int *ld=LOGICAL(data);
				for(i=0;i<xlength(data);i++) {
					bdata[i]=ld[i]==NA_LOGICAL ? 255 : ld[i];
					}
				offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, xlength(data), bdata, *moffset);
				free(bdata);
				break;
			case INTSXP:
				bdata=calloc(xlength(data), 1);
				if(bdata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				int *id=INTEGER(data);
				for(i=0;i<xlength(data);i++)bdata[i]=id[i];
				offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, xlength(data), bdata, *moffset);
				free(bdata);
				break;
			case REALSXP:
				bdata=calloc(xlength(data), 1);
				if(bdata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				double *fd=REAL(data);
				for(i=0;i<xlength(data);i++)bdata[i]=fd[i];
				offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, xlength(data), bdata, *moffset);
				free(bdata);
				break;
			case STRSXP:
				if(xlength(data)!=1) {
					error("Can only convert a single string to UINT8");
					return(R_NilValue);
					}
				SEXP sch=STRING_ELT(data, 0);
				if(sch==NA_STRING) {
					offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, MVL_NA_STRING_LENGTH, MVL_NA_STRING, *moffset);
					} else {
					const char *bd=CHAR(sch);
					offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, strlen(bd), bd, *moffset);
					}
				break;
			default:
				error("Cannot convert R type %d to UINT8", TYPEOF(data));
				return(R_NilValue);
				break;
			}
		break;
	case LIBMVL_VECTOR_INT32:
		offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT32, xlength(data), INTEGER(data), *moffset);
		break;
	case LIBMVL_VECTOR_INT64:
		switch(TYPEOF(data)) {
			case RAWSXP:
//				Rprintf("Writing INT64 from RAW input length %lld output length %lld\n", xlength(data), xlength(data)/8);
				offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, xlength(data)/8, RAW(data), *moffset);
				break;
			case REALSXP:
				idata=calloc(xlength(data), sizeof(*idata));
				if(idata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				pd=REAL(data);
				for(i=0;i<xlength(data);i++)
					idata[i]=pd[i];
				
				offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, xlength(data), idata, *moffset);
				free(idata);
				break;
			case INTSXP:
				idata=calloc(xlength(data), sizeof(*idata));
				if(idata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				pi=INTEGER(data);
				for(i=0;i<xlength(data);i++)
					idata[i]=pi[i];
				
				offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, xlength(data), idata, *moffset);
				free(idata);
				break;
			default:
				error("can only write raw, double and integer to INT64");
				return(R_NilValue);
			}
		break;
	case LIBMVL_VECTOR_FLOAT:
		fdata=calloc(xlength(data), sizeof(*fdata));
		if(fdata==NULL) {
			error("Out of memory");
			return(R_NilValue);
			}
		pd=REAL(data);
		for(i=0;i<xlength(data);i++)
			fdata[i]=pd[i];
		offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_FLOAT, xlength(data), fdata, *moffset);
		free(fdata);
		break;
	case LIBMVL_VECTOR_DOUBLE:
		offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_DOUBLE, xlength(data), REAL(data), *moffset);
		break;
	case LIBMVL_VECTOR_OFFSET64:
		offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_OFFSET64, xlength(data), REAL(data), *moffset);
		break;
#if 0
	case 10000:
		strvec=calloc(xlength(data), sizeof(*strvec));
		if(strvec==NULL) {
			error("Out of memory");
			return(R_NilValue);
			}
		for(i=0;i<xlength(data);i++) {
			ch=CHAR(STRING_ELT(data, i));
			strvec[i]=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_CSTRING, strlen(ch), ch, 0);
			}
		offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_OFFSET64, xlength(data), strvec, *moffset);
		free(strvec);
		break;
#else
	case 10000:
		strvec2=calloc(xlength(data), sizeof(*strvec2));
		strlen2=calloc(xlength(data), sizeof(*strlen2));
		if(strvec2==NULL || strlen2==NULL) {
			error("Out of memory");
			free(strvec2);
			free(strlen2);
			return(R_NilValue);
			}
		for(i=0;i<xlength(data);i++) {
			SEXP ch=STRING_ELT(data, i);
			if(ch==NA_STRING) {
				strvec2[i]=MVL_NA_STRING;
				strlen2[i]=MVL_NA_STRING_LENGTH;
				} else {
				strvec2[i]=CHAR(ch);
				strlen2[i]=xlength(ch);
				}
			}
		offset=mvl_write_packed_list(libraries[idx].ctx, xlength(data), strlen2, (unsigned char **)strvec2, *moffset);
		free(strvec2);
		free(strlen2);
		break;
#endif
	case 10001:
		if(length(data)!=1) {
			error("data has to be length 1 string vector");
			return(R_NilValue);
			}
		SEXP sch=STRING_ELT(data, 0);
		if(sch==NA_STRING) 
			offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_CSTRING, MVL_NA_STRING_LENGTH, MVL_NA_STRING, *moffset);
			else {
			ch=CHAR(sch);
			offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_CSTRING, strlen(ch), ch, *moffset);
			}
		break;
		
	default:
		error("write_vector: unknown type %d", type);
		return(R_NilValue);
	}
ans=PROTECT(allocVector(REALSXP, 1));
REAL(ans)[0]=*doffset;

class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(2);
return(ans);
}

SEXP start_write_vector(SEXP idx0, SEXP type0, SEXP expected_length0, SEXP data, SEXP metadata_offset)
{
long i;
int idx, type;
double dmoffset;
LIBMVL_OFFSET64 *moffset=(LIBMVL_OFFSET64 *)&dmoffset;

LIBMVL_OFFSET64 offset, expected_length;
double *doffset=(double *)&offset;
const char *ch, **strvec2;
LIBMVL_OFFSET64 *strvec;
long long *idata;
float *fdata;
unsigned char *bdata;

double *pd;
int *pi;

SEXP ans, class;

if(length(idx0)!=1) {
	error("write_vector first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].f==NULL) {
	error("library not open for writing");
	return(R_NilValue);
	}
if(mvl_get_error(libraries[idx].ctx)!=0) {
	error("library has error status %d: %s", mvl_get_error(libraries[idx].ctx), mvl_strerror(libraries[idx].ctx));
	return(R_NilValue);
	}

if(length(type0)!=1) {
	error("write_vector second argument must be a single integer");
	return(R_NilValue);
	}
type=INTEGER(type0)[0];

if(length(expected_length0)!=1) {
	error("third parameter parameter (expected_length0) must be a single real or integer");
	return(R_NilValue);
	}
	
switch(TYPEOF(expected_length0)) {
	case REALSXP:
		expected_length=REAL(expected_length0)[0];
		break;
	case INTSXP:
		expected_length=INTEGER(expected_length0)[0];
		break;
	case NILSXP:
		expected_length=xlength(data);
		break;
	default:
		error("Unknown expected offset R type %d", TYPEOF(expected_length0));
		return(R_NilValue);
	}


libraries[idx].modified=1;

if(length(metadata_offset)<1) {
	*moffset=0;
	} else {
	dmoffset=REAL(metadata_offset)[0];
	}
switch(type) {
	case LIBMVL_VECTOR_UINT8:
		switch(TYPEOF(data)) {
			case RAWSXP:
				offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, expected_length, xlength(data), RAW(data), *moffset);
				break;
			case INTSXP:
				bdata=calloc(xlength(data), 1);
				if(bdata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				int *id=INTEGER(data);
				for(i=0;i<xlength(data);i++)bdata[i]=id[i];
				offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, expected_length, xlength(data), bdata, *moffset);
				free(bdata);
				break;
			case LGLSXP:
				bdata=calloc(xlength(data), 1);
				if(bdata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				int *ld=LOGICAL(data);
				for(i=0;i<xlength(data);i++) {
					bdata[i]=ld[i]==NA_LOGICAL ? 255 : ld[i];
					}
				offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, expected_length, xlength(data), bdata, *moffset);
				free(bdata);
				break;
			case REALSXP:
				bdata=calloc(xlength(data), 1);
				if(bdata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				double *fd=REAL(data);
				for(i=0;i<xlength(data);i++)bdata[i]=fd[i];
				offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, expected_length, xlength(data), bdata, *moffset);
				free(bdata);
				break;
			case STRSXP:
				if(xlength(data)!=1) {
					error("Can only convert a single string to UINT8");
					return(R_NilValue);
					}
				SEXP ch=STRING_ELT(data, 0);
				if(ch==NA_STRING)
					offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, expected_length, MVL_NA_STRING_LENGTH, MVL_NA_STRING, *moffset);
					else {
					const char *bd=CHAR(ch);
					offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, expected_length, strlen(bd), bd, *moffset);
					}
				break;
			default:
				error("Cannot convert R type %d to UINT8", TYPEOF(data));
				return(R_NilValue);
				break;
			}
		break;
	case LIBMVL_VECTOR_INT32:
		offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT32, expected_length, xlength(data), INTEGER(data), *moffset);
		break;
	case LIBMVL_VECTOR_INT64:
		switch(TYPEOF(data)) {
			case RAWSXP:
//				Rprintf("Writing INT64 from RAW input length %lld output length %lld\n", xlength(data), xlength(data)/8);
				offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, expected_length, xlength(data)/8, RAW(data), *moffset);
				break;
			case REALSXP:
				idata=calloc(xlength(data), sizeof(*idata));
				if(idata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				pd=REAL(data);
				for(i=0;i<xlength(data);i++)
					idata[i]=pd[i];
				
				offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, expected_length, xlength(data), idata, *moffset);
				free(idata);
				break;
			case INTSXP:
				idata=calloc(xlength(data), sizeof(*idata));
				if(idata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				pi=INTEGER(data);
				for(i=0;i<xlength(data);i++)
					idata[i]=pi[i];
				
				offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, expected_length, xlength(data), idata, *moffset);
				free(idata);
				break;
			default:
				error("can only write raw, double and integer to INT64");
				return(R_NilValue);
			}
		break;
	case LIBMVL_VECTOR_FLOAT:
		fdata=calloc(xlength(data), sizeof(*fdata));
		if(fdata==NULL) {
			error("Out of memory");
			return(R_NilValue);
			}
		pd=REAL(data);
		for(i=0;i<xlength(data);i++)
			fdata[i]=pd[i];
		offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_FLOAT, expected_length, xlength(data), fdata, *moffset);
		free(fdata);
		break;
	case LIBMVL_VECTOR_DOUBLE:
		offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_DOUBLE, expected_length, xlength(data), REAL(data), *moffset);
		break;
	case LIBMVL_VECTOR_OFFSET64:
		offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_OFFSET64, expected_length, xlength(data), REAL(data), *moffset);
		break;
#if 0
	case 10000:
		strvec=calloc(xlength(data), sizeof(*strvec));
		if(strvec==NULL) {
			error("Out of memory");
			return(R_NilValue);
			}
		for(i=0;i<xlength(data);i++) {
			ch=CHAR(STRING_ELT(data, i));
			strvec[i]=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_CSTRING, strlen(ch), ch, 0);
			}
		offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_OFFSET64, xlength(data), strvec, *moffset);
		free(strvec);
		break;
#else
	case 10000:
		error("Rewriting packed lists is not supported");
		return(R_NilValue);
// 		strvec2=calloc(xlength(data), sizeof(*strvec2));
// 		if(strvec2==NULL) {
// 			error("Out of memory");
// 			return(R_NilValue);
// 			}
// 		for(i=0;i<xlength(data);i++) {
// 			strvec2[i]=CHAR(STRING_ELT(data, i));
// 			}
// 		offset=mvl_write_packed_list(libraries[idx].ctx, expected_length, xlength(data), NULL, (char **)strvec2, *moffset);
// 		free(strvec2);
// 		break;
#endif
	case 10001:
		if(length(data)!=1) {
			error("data has to be length 1 string vector");
			return(R_NilValue);
			}
		SEXP sch=STRING_ELT(data, 0);
		if(sch==NA_STRING)
			offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_CSTRING, expected_length, MVL_NA_STRING_LENGTH, MVL_NA_STRING, *moffset);
			else {
			ch=CHAR(sch);
			offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_CSTRING, expected_length, strlen(ch), ch, *moffset);
			}
		break;
		
	default:
		error("write_vector: unknown type %d", type);
		return(R_NilValue);
	}
ans=PROTECT(allocVector(REALSXP, 1));
REAL(ans)[0]=*doffset;

class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(2);
return(ans);
}

SEXP rewrite_vector(SEXP vec0, SEXP offset, SEXP data)
{
long i;
int idx, type;

LIBMVL_OFFSET64 data_offset, Nelements, chunk_offset;
const char *ch, **strvec2;
LIBMVL_OFFSET64 *strvec;
long long *idata;
float *fdata;
unsigned char *bdata;
double *ddata;

double *pd;
int *pi;
LIBMVL_VECTOR *vec;

decode_mvl_object(vec0, &idx, &data_offset);
vec=get_mvl_vector(idx, data_offset);

if(vec==NULL) {
	error("Invalid MVL object");
	return(R_NilValue);
	}

if(libraries[idx].f==NULL) {
	error("library not open for writing");
	return(R_NilValue);
	}
if(mvl_get_error(libraries[idx].ctx)!=0) {
	error("library has error status %d: %s", mvl_get_error(libraries[idx].ctx), mvl_strerror(libraries[idx].ctx));
	return(R_NilValue);
	}
	
if(length(offset)!=1) {
	error("second parameter (offset) must be a single real or integer");
	return(R_NilValue);
	}

switch(TYPEOF(offset)) {
	case REALSXP:
		chunk_offset=REAL(offset)[0]-1;
		break;
	case INTSXP:
		chunk_offset=INTEGER(offset)[0]-1;
		break;
	default:
		error("Unknown offset R type %d", TYPEOF(offset));
		return(R_NilValue);
	}

type=mvl_vector_type(vec);


switch(type) {
	case LIBMVL_PACKED_LIST64:
		Nelements=mvl_vector_length(vec)-1;
		break;
	default:
		Nelements=mvl_vector_length(vec);
		break;		
	}
	
if(chunk_offset+xlength(data)>Nelements) {
	error("Attempt to write past vector end");
	return(R_NilValue);
	}
if(chunk_offset<0) {
	error("Attempt to write before index 1");
	return(R_NilValue);
	}
	

libraries[idx].modified=1;

switch(type) {
	case LIBMVL_VECTOR_UINT8:
		switch(TYPEOF(data)) {
			case RAWSXP:
				mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data), RAW(data));
				break;
			case INTSXP:
				bdata=calloc(xlength(data), 1);
				if(bdata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				int *id=INTEGER(data);
				for(i=0;i<xlength(data);i++)bdata[i]=id[i];
				mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data), bdata);
				free(bdata);
				break;
			case LGLSXP:
				bdata=calloc(xlength(data), 1);
				if(bdata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				int *ld=LOGICAL(data);
				for(i=0;i<xlength(data);i++) {
					bdata[i]=ld[i]==NA_LOGICAL ? 255 : ld[i];
					}
				mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data), bdata);
				free(bdata);
				break;
			case REALSXP:
				bdata=calloc(xlength(data), 1);
				if(bdata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				double *fd=REAL(data);
				for(i=0;i<xlength(data);i++)bdata[i]=fd[i];
				mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data), bdata);
				free(bdata);
				break;
			case STRSXP:
				if(xlength(data)!=1) {
					error("Can only convert a single string to UINT8");
					return(R_NilValue);
					}
				SEXP ch=STRING_ELT(data, 0);
				if(ch==NA_STRING)
					mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, MVL_NA_STRING_LENGTH, MVL_NA_STRING);
					else {
					const char *bd=CHAR(ch);
					mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, strlen(bd), bd);
					}
				break;
			default:
				error("Cannot convert R type %d to UINT8", TYPEOF(data));
				return(R_NilValue);
				break;
			}
		break;
	case LIBMVL_VECTOR_INT32:
		mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data), INTEGER(data));
		break;
	case LIBMVL_VECTOR_INT64:
		switch(TYPEOF(data)) {
			case RAWSXP:
//				Rprintf("Writing INT64 from RAW input length %lld output length %lld\n", xlength(data), xlength(data)/8);
				mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data)/8, RAW(data));
				break;
			case REALSXP:
				idata=calloc(xlength(data), sizeof(*idata));
				if(idata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				pd=REAL(data);
				for(i=0;i<xlength(data);i++)
					idata[i]=pd[i];
				
				mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data), idata);
				free(idata);
				break;
			case INTSXP:
				idata=calloc(xlength(data), sizeof(*idata));
				if(idata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				pi=INTEGER(data);
				for(i=0;i<xlength(data);i++)
					idata[i]=pi[i];
				
				mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data), idata);
				free(idata);
				break;
			default:
				error("can only write raw, double and integer to INT64");
				return(R_NilValue);
			}
		break;
	case LIBMVL_VECTOR_FLOAT:
		fdata=calloc(xlength(data), sizeof(*fdata));
		if(fdata==NULL) {
			error("Out of memory");
			return(R_NilValue);
			}
		switch(TYPEOF(data)) {
			case REALSXP:
				pd=REAL(data);
				for(i=0;i<xlength(data);i++)
					fdata[i]=pd[i];
				break;
			case INTSXP:
				pi=INTEGER(data);
				for(i=0;i<xlength(data);i++)
					fdata[i]=pi[i];
				break;
			default:
				error("Cannot convert R type %d to floating point", TYPEOF(data));
				return(R_NilValue);
			}
				
		mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data), fdata);
		free(fdata);
		break;
	case LIBMVL_VECTOR_DOUBLE:
		switch(TYPEOF(data)) {
			case REALSXP:
				mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data), REAL(data));
				break;
			case INTSXP: {
				ddata=calloc(xlength(data), sizeof(*ddata));
				if(ddata==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				pi=INTEGER(data);
				for(i=0;i<xlength(data);i++)
					ddata[i]=pi[i];
				mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data), ddata);
				free(ddata);
				break;
				}
			default:
				error("Cannot convert R type %d to floating point", TYPEOF(data));
				return(R_NilValue);
			}
		break;
	case LIBMVL_VECTOR_OFFSET64:
		mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, xlength(data), REAL(data));
		break;
#if 0
	case 10000:
		strvec=calloc(xlength(data), sizeof(*strvec));
		if(strvec==NULL) {
			error("Out of memory");
			return(R_NilValue);
			}
		for(i=0;i<xlength(data);i++) {
			ch=CHAR(STRING_ELT(data, i));
			strvec[i]=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_CSTRING, strlen(ch), ch, 0);
			}
		offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_OFFSET64, xlength(data), strvec, *moffset);
		free(strvec);
		break;
#else
		/* Rewriting packed lists is not supported */
// 	case 10000:
// 		strvec2=calloc(xlength(data), sizeof(*strvec2));
// 		if(strvec2==NULL) {
// 			error("Out of memory");
// 			return(R_NilValue);
// 			}
// 		for(i=0;i<xlength(data);i++) {
// 			strvec2[i]=CHAR(STRING_ELT(data, i));
// 			}
// 		mvl_rewrite_packed_list(libraries[idx].ctx, data_offset, chunk_offset, xlength(data), NULL, (char **)strvec2);
// 		free(strvec2);
// 		break;
#endif
	case LIBMVL_VECTOR_CSTRING:
		if(length(data)!=1) {
			error("data has to be length 1 string vector");
			return(R_NilValue);
			}
		SEXP sch=STRING_ELT(data, 0);
		if(sch==NA_STRING)
			mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, MVL_NA_STRING_LENGTH, MVL_NA_STRING);
			else {
			ch=CHAR(sch);
			mvl_rewrite_vector(libraries[idx].ctx, type, data_offset, chunk_offset, strlen(ch), ch);
			}
		break;
		
	default:
		error("write_vector: unsupported type %d", type);
		return(R_NilValue);
	}
return(R_NilValue);
}

SEXP fused_write_vector(SEXP idx0, SEXP type0, SEXP data_list, SEXP metadata_offset)
{
long i, j, k, m;
int idx, data_idx, type;
double dmoffset;
LIBMVL_OFFSET64 *moffset=(LIBMVL_OFFSET64 *)&dmoffset;

LIBMVL_OFFSET64 offset, char_offset, data_offset, vec_idx, char_idx;
double *doffset=(double *)&offset;
const char *ch;
LIBMVL_OFFSET64 *strvec;
long long *idata;
float *fdata;
SEXP data, sch;

double *pd;
int *pi;
unsigned char *pc;
long total_length, char_total_length;

LIBMVL_VECTOR *vec;

SEXP ans, class;

if(length(idx0)!=1) {
	error("fused_write_vector first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].f==NULL) {
	error("library not open for writing");
	return(R_NilValue);
	}
if(mvl_get_error(libraries[idx].ctx)!=0) {
	error("library has error status %d: %s", mvl_get_error(libraries[idx].ctx), mvl_strerror(libraries[idx].ctx));
	return(R_NilValue);
	}
	
if(TYPEOF(data_list)!=VECSXP) {
	error("fused_write_vector third argument must be a list of data to write out");
	return(R_NilValue);
	}

if(length(type0)!=1) {
	error("fused_write_vector second argument must be a single integer");
	return(R_NilValue);
	}
type=INTEGER(type0)[0];

libraries[idx].modified=1;

if(length(metadata_offset)<1) {
	*moffset=0;
	} else {
	dmoffset=REAL(metadata_offset)[0];
	}

total_length=0;
char_total_length=0;
for(k=0;k<xlength(data_list);k++) {
	data=PROTECT(VECTOR_ELT(data_list, k));
	
	if(TYPEOF(data)==STRSXP) {
		total_length+=xlength(data);
		for(i=0;i<xlength(data);i++) {
			sch=STRING_ELT(data, i);
			if(sch==NA_STRING)
				char_total_length+=MVL_NA_STRING_LENGTH;
				else
				char_total_length+=strlen(CHAR(sch));
			}
		UNPROTECT(1);
		continue;
		}
	
	if(TYPEOF(data)!=VECSXP) {
		total_length+=xlength(data);
		UNPROTECT(1);
		continue;
		}
	decode_mvl_object(data, &data_idx, &data_offset);
	UNPROTECT(1);
	vec=get_mvl_vector(data_idx, data_offset);
	
	if(vec==NULL) {
		error("Invalid MVL object in data list");
		return(R_NilValue);
		}
	if(mvl_vector_type(vec)!=type) {
		switch(type) {
			case 10000:
				if(mvl_vector_type(vec)!=LIBMVL_PACKED_LIST64) {
					error("Cannot convert MVL object to character");
					return(R_NilValue);
					}
				if(mvl_vector_length(vec)<1) {
					error("Invalid MVL packed list");
					return(R_NilValue);
					}
				//Rprintf("STRVEC length %ld\n", mvl_vector_length(vec));
				total_length+=mvl_vector_length(vec)-1;
				char_total_length+=mvl_vector_data_offset(vec)[mvl_vector_length(vec)-1]-mvl_vector_data_offset(vec)[0];
				break;
			default:
				error("Internal conversion between types of MVL objects not supported yet");
				return(R_NilValue);
			}
		} else {
		total_length+=mvl_vector_length(vec);
		}
	}

if(char_total_length>0) {
	char_offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, char_total_length, 0, NULL, LIBMVL_NO_METADATA);
	char_offset+=sizeof(LIBMVL_VECTOR_HEADER);
	total_length++;
	offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_PACKED_LIST64, total_length, 1, &char_offset, *moffset);
	//Rprintf("packed_list: %ld %ld %ld %ld\n", total_length, char_total_length, offset, char_offset);
	vec_idx=1;
	char_idx=0;
	} else {
	offset=mvl_start_write_vector(libraries[idx].ctx, type, total_length, 0, NULL, *moffset);
	vec_idx=0;
	char_idx=0;
	}

for(k=0;k<xlength(data_list);k++) {
	data=PROTECT(VECTOR_ELT(data_list, k));

	if(TYPEOF(data)==VECSXP) {
		decode_mvl_object(data, &data_idx, &data_offset);
		vec=get_mvl_vector(data_idx, data_offset);
		
		if(mvl_vector_type(vec)==type) {
			mvl_rewrite_vector(libraries[idx].ctx, type, offset, vec_idx, mvl_vector_length(vec), (mvl_vector_data_uint8(vec)));
			vec_idx+=mvl_vector_length(vec);
			UNPROTECT(1);
			continue;
			}

		switch(type) {
			case LIBMVL_VECTOR_INT32:
				break;
			case LIBMVL_VECTOR_INT64:
				break;
			case LIBMVL_VECTOR_FLOAT:
				break;
			case LIBMVL_VECTOR_DOUBLE:
				break;
			case 10000:
				pc=sizeof(LIBMVL_VECTOR_HEADER)+(unsigned char *)get_mvl_vector(data_idx, mvl_vector_data_offset(vec)[0]-sizeof(LIBMVL_VECTOR_HEADER));
				if(pc==NULL) {
					error("Invalid packed list");
					return(R_NilValue);
					}
				mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, char_offset-sizeof(LIBMVL_VECTOR_HEADER), char_idx, mvl_vector_data_offset(vec)[mvl_vector_length(vec)-1]-mvl_vector_data_offset(vec)[0], pc);
				//Rprintf("s# %ld %ld\n", vec_idx, char_idx);

				#define REWRITE_BUF_SIZE (1024*1024)
				
				strvec=calloc(REWRITE_BUF_SIZE, sizeof(*strvec));
				if(strvec==NULL) {
					error("Out of memory");
					return(R_NilValue);
					}
				
				for(j=1;j<mvl_vector_length(vec);j+=REWRITE_BUF_SIZE) {
					for(i=0;i+j<mvl_vector_length(vec) && i<REWRITE_BUF_SIZE;i++)
						strvec[i]=char_offset+char_idx+(mvl_vector_data_offset(vec)[j+i]-mvl_vector_data_offset(vec)[0]);
					i=j+REWRITE_BUF_SIZE>=mvl_vector_length(vec) ? mvl_vector_length(vec)-j : REWRITE_BUF_SIZE;
					mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_PACKED_LIST64, offset, vec_idx, i, strvec);
					vec_idx+=i;
					}
					
				
				/* TODO: all offsets need to shift by char_idx-mvl_vector_data_offset(vec)[0] */
				
				char_idx+=mvl_vector_data_offset(vec)[mvl_vector_length(vec)-1]-mvl_vector_data_offset(vec)[0];
				//Rprintf("s#2 %ld %ld\n", vec_idx, char_idx);
				
				free(strvec);
				#undef REWRITE_BUF_SIZE
				break;
			default:
				break;
			}
		/* TODO */
		UNPROTECT(1);
		continue;
		}
	
	switch(type) {
		case LIBMVL_VECTOR_UINT8:
			switch(TYPEOF(data)) {
				case RAWSXP:
					mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, offset, vec_idx, xlength(data), RAW(data));
					break;
				case LGLSXP:
					pc=malloc(xlength(data));
					if(pc==NULL) {
						error("Out of memory");
						return(R_NilValue);
						}
					int *ld=LOGICAL(data);
					for(i=0;i<xlength(data);i++) {
						pc[i]=ld[i]==NA_LOGICAL ? 255 : ld[i];
						}
					mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, offset, vec_idx, xlength(data), pc);
					free(pc);
					break;
				default:
					error("Cannot convert data with type=%d to raw vector", TYPEOF(data));
					return(R_NilValue);
					break;
				}
			vec_idx+=xlength(data);
			break;
		case LIBMVL_VECTOR_INT32:
			mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT32, offset, vec_idx, xlength(data), INTEGER(data));
			vec_idx+=xlength(data);
			break;
		case LIBMVL_VECTOR_INT64:
			switch(TYPEOF(data)) {
				case RAWSXP:
	//				Rprintf("Writing INT64 from RAW input length %lld output length %lld\n", xlength(data), xlength(data)/8);
					mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, offset, vec_idx, xlength(data)/8, RAW(data));
					vec_idx+=xlength(data)/8;
					break;
				case REALSXP:
					idata=calloc(xlength(data), sizeof(*idata));
					if(idata==NULL) {
						error("Out of memory");
						return(R_NilValue);
						}
					pd=REAL(data);
					for(i=0;i<xlength(data);i++)
						idata[i]=pd[i];
					
					mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, offset, vec_idx, xlength(data), idata);
					vec_idx+=xlength(data);
					free(idata);
					break;
				case INTSXP:
					idata=calloc(xlength(data), sizeof(*idata));
					if(idata==NULL) {
						error("Out of memory");
						return(R_NilValue);
						}
					pi=INTEGER(data);
					for(i=0;i<xlength(data);i++)
						idata[i]=pi[i];
					
					mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, offset, vec_idx, xlength(data), idata);
					vec_idx+=xlength(data);
					free(idata);
					break;
				default:
					error("can only write raw, double and integer to INT64");
					UNPROTECT(1);
					return(R_NilValue);
				}
			break;
		case LIBMVL_VECTOR_FLOAT:
			fdata=calloc(xlength(data), sizeof(*fdata));
			if(fdata==NULL) {
				error("Out of memory");
				UNPROTECT(1);
				return(R_NilValue);
				}
			pd=REAL(data);
			for(i=0;i<xlength(data);i++)
				fdata[i]=pd[i];
			mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_FLOAT, offset, vec_idx, xlength(data), fdata);
			vec_idx+=xlength(data);
			free(fdata);
			break;
		case LIBMVL_VECTOR_DOUBLE:
			mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_DOUBLE, offset, vec_idx, xlength(data), REAL(data));
			vec_idx+=xlength(data);
			break;
		case LIBMVL_VECTOR_OFFSET64:
			mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_OFFSET64, offset, vec_idx, xlength(data), REAL(data));
			vec_idx+=xlength(data);
			break;
		case 10000:
			#define REWRITE_BUF_SIZE (1024*1024)
			
			strvec=calloc(REWRITE_BUF_SIZE, sizeof(*strvec));
			if(strvec==NULL) {
				error("Out of memory");
				return(R_NilValue);
				}
			
			/* TODO: it would be nice to bounce strings against internal buffer to reduce frequency of rewrite() calls */
			for(j=0;j<xlength(data);j+=REWRITE_BUF_SIZE) {
				for(i=0;i+j<xlength(data) && i<REWRITE_BUF_SIZE;i++) {
					sch=STRING_ELT(data, i+j);
					if(sch==NA_STRING) {
						ch=MVL_NA_STRING;
						m=MVL_NA_STRING_LENGTH;
						mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, char_offset-sizeof(LIBMVL_VECTOR_HEADER), char_idx, m, ch);
						//Rprintf("str %s %d %ld\n", ch, m, char_idx);
						strvec[i]=char_offset+char_idx+m;
						char_idx+=m;
						} else {
						ch=CHAR(sch);
						m=strlen(ch);
						mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_UINT8, char_offset-sizeof(LIBMVL_VECTOR_HEADER), char_idx, m, ch);
						//Rprintf("str %s %d %ld\n", ch, m, char_idx);
						strvec[i]=char_offset+char_idx+m;
						char_idx+=m;
						}
					}
				i=j+REWRITE_BUF_SIZE>=xlength(data) ? xlength(data)-j : REWRITE_BUF_SIZE;
				mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_PACKED_LIST64, offset, vec_idx, i, strvec);
				vec_idx+=i;
				}
			#undef REWRITE_BUF_SIZE
			free(strvec);				
			break;
// 		case 10001:
// 			if(length(data)!=1) {
// 				error("data has to be length 1 string vector");
// 				return(R_NilValue);
// 				}
// 			ch=CHAR(STRING_ELT(data, 0));
// 			offset=mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_CSTRING, strlen(ch), ch, *moffset);
// 			break;
			
		default:
			error("write_vector: unknown type %d", type);
			UNPROTECT(1);
			return(R_NilValue);
		}
	UNPROTECT(1);
	}
	
if(char_idx>char_total_length) {
	error("Internal error: char_idx>char_total_length");
	}
	
ans=PROTECT(allocVector(REALSXP, 1));
REAL(ans)[0]=*doffset;

class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(2);
return(ans);
}

SEXP indexed_copy_vector(SEXP idx0, SEXP mvl_object, SEXP indices, SEXP metadata_offset)
{
LIBMVL_OFFSET64 *v_idx;
LIBMVL_OFFSET64 offset, data_offset, N_idx;
double *doffset=(double *)&offset;
int idx, data_idx;

double dmoffset;
LIBMVL_OFFSET64 *moffset=(LIBMVL_OFFSET64 *)&dmoffset;

LIBMVL_VECTOR *vec;

SEXP ans, class;

if(length(idx0)!=1) {
	error("fused_write_vector first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].f==NULL) {
	error("library not open for writing");
	return(R_NilValue);
	}
if(mvl_get_error(libraries[idx].ctx)!=0) {
	error("library has error status %d: %s", mvl_get_error(libraries[idx].ctx), mvl_strerror(libraries[idx].ctx));
	return(R_NilValue);
	}
		
if(length(metadata_offset)<1) {
	*moffset=0;
	} else {
	dmoffset=REAL(metadata_offset)[0];
	}
	
if(TYPEOF(mvl_object)!=VECSXP) {
	error("Not a valid MVL object");
	return(R_NilValue);	
	}
	
decode_mvl_object(mvl_object, &data_idx, &data_offset);
vec=get_mvl_vector(data_idx, data_offset);

if(vec==NULL) {
	error("Not a valid MVL object (2)");
	return(R_NilValue);	
	}

	
if(get_indices(indices, vec, &N_idx, &v_idx)) {
	error("Invalid indices");
	return(R_NilValue);
	}
	
	
libraries[idx].modified=1;

offset=mvl_indexed_copy_vector(libraries[idx].ctx, N_idx, v_idx, vec, libraries[data_idx].data, libraries[data_idx].length, *moffset, 1024*1024*16);
	
free(v_idx);

ans=PROTECT(allocVector(REALSXP, 1));
REAL(ans)[0]=*doffset;

class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(2);

return(ans);
}

SEXP order_vectors(SEXP data_list, SEXP indices, SEXP s_sort_function)
{
int data_idx, sort_function;
LIBMVL_OFFSET64 data_offset;

double *pd;
int *pi, err;

void **vec_data;
LIBMVL_VECTOR **vectors, *vec;
LIBMVL_OFFSET64 *v_idx, *vec_data_size;
LIBMVL_OFFSET64 N;
LIBMVL_OFFSET64 vector_length;

SEXP ans;
	
if(TYPEOF(data_list)!=VECSXP) {
	error("order_vectors first argument must be a list of data to sort");
	return(R_NilValue);
	}

if(xlength(data_list)<1) {
	return(indices);
	}
	
if(TYPEOF(indices)!=NILSXP && xlength(indices)<1) {
	return(indices);
	}
	
if(TYPEOF(s_sort_function)!=INTSXP || xlength(s_sort_function)!=1) {
	error("Invalid sort function");
	return(R_NilValue);
	}

sort_function=INTEGER(s_sort_function)[0];
	
vec_data=calloc(xlength(data_list), sizeof(*vec_data));
vec_data_size=calloc(xlength(data_list), sizeof(*vec_data_size));
vectors=calloc(xlength(data_list), sizeof(*vectors));
if(vec_data==NULL || vectors==NULL || vec_data_size==NULL) {
	error("Not enough memory");
	free(vec_data);
	free(vec_data_size);
	free(vectors);
	return(R_NilValue);
	}

vector_length=0;
for(LIBMVL_OFFSET64 k=0;k<xlength(data_list);k++) {
	decode_mvl_object(PROTECT(VECTOR_ELT(data_list, k)), &data_idx, &data_offset);
	UNPROTECT(1);
	vectors[k]=get_mvl_vector(data_idx, data_offset);
	
	if(vectors[k]==NULL) {
		error("Invalid MVL object in data list");
		free(vec_data);
		free(vec_data_size);
		free(vectors);
		return(R_NilValue);
		}
		
	if(k==0) { 
		vector_length=mvl_vector_nentries(vectors[k]);
		} else {
		if(vector_length!=mvl_vector_nentries(vectors[k])) {
			error("Mismatched length %llu (expected %llu) for vector %llu", mvl_vector_length(vectors[k]), vector_length, k);
			free(vec_data);
			free(vec_data_size);
			free(vectors);
			return(R_NilValue);
			}
		}
		
	vec_data[k]=libraries[data_idx].data;
	vec_data_size[k]=libraries[data_idx].length;
	}
	
switch(TYPEOF(indices)) {
	case VECSXP:
		decode_mvl_object(indices, &data_idx, &data_offset);
		vec=get_mvl_vector(data_idx, data_offset);
		
		if(vec==NULL) {
			error("Invalid MVL object or R vector passed as indices");
			free(vec_data);
			free(vec_data_size);
			free(vectors);
			return(R_NilValue);
			}
		N=mvl_vector_length(vec);
		v_idx=calloc(N, sizeof(*v_idx));
		if(v_idx==NULL) {
			error("Not enough memory");
			free(vec_data);
			free(vec_data_size);
			free(vectors);
			return(R_NilValue);
			}
		
		switch(mvl_vector_type(vec)) {
			case LIBMVL_VECTOR_OFFSET64:
				memcpy(v_idx, mvl_vector_data_offset(vec), N*sizeof(*v_idx));
				break;
			case LIBMVL_VECTOR_INT32:
				for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=mvl_vector_data_int32(vec)[m];
				break;
			case LIBMVL_VECTOR_INT64:
				for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=mvl_vector_data_int64(vec)[m];
				break;
			case LIBMVL_VECTOR_DOUBLE:
				for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=mvl_vector_data_double(vec)[m];
				break;
			case LIBMVL_VECTOR_FLOAT:
				for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=mvl_vector_data_float(vec)[m];
				break;
			default:
				error("Cannot interpret MVL object as indices");
				free(vec_data);
				free(vec_data_size);
				free(vectors);
				free(v_idx);
				return(R_NilValue);
				break;
			}
		break;
	case REALSXP:
		N=xlength(indices);
		v_idx=calloc(N, sizeof(*v_idx));
		if(v_idx==NULL) {
			error("Not enough memory");
			free(vec_data);
			free(vec_data_size);
			free(vectors);
			return(R_NilValue);
			}
		pd=REAL(indices);
		for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=pd[m]-1;
		break;
	case INTSXP:
		N=xlength(indices);
		v_idx=calloc(N, sizeof(*v_idx));
		if(v_idx==NULL) {
			error("Not enough memory");
			free(vec_data);
			free(vec_data_size);
			free(vectors);
			return(R_NilValue);
			}
		pi=INTEGER(indices);
		for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=pi[m]-1;
		break;
	case NILSXP:
		N=mvl_vector_length(vectors[0]);
		if(mvl_vector_type(vectors[0])==LIBMVL_PACKED_LIST64)N--;
		v_idx=calloc(N, sizeof(*v_idx));
		if(v_idx==NULL) {
			error("Not enough memory");
			free(vec_data);
			free(vec_data_size);
			free(vectors);
			return(R_NilValue);
			}
		for(LIBMVL_OFFSET64 m=0;m<N;m++)v_idx[m]=m;
		break;
	default:
		error("Cannot interpret R object as index");
		free(vec_data);
		free(vec_data_size);
		free(vectors);
		return(R_NilValue);		
	}
	
for(LIBMVL_OFFSET64 m=0;m<N;m++) {
	if(v_idx[m]>=vector_length) {
		error("Index %llu larger than vector length (%llu vs %llu)", m+1, v_idx[m]+1, vector_length);
		free(v_idx);
		free(vec_data);
		free(vec_data_size);
		free(vectors);
		return(R_NilValue);				
		}
	}
	
for(LIBMVL_OFFSET64 k=0;k<xlength(data_list);k++) {
	if(mvl_vector_type(vectors[k])!=LIBMVL_PACKED_LIST64)continue;
	
	
	for(LIBMVL_OFFSET64 m=0;m<N;m++) {
		if(mvl_packed_list_validate_entry(vectors[k], vec_data[k], vec_data_size[k], v_idx[m])!=0) {
			error("Invalid entry %llu in string vector %llu", m+1, k+1);
			free(v_idx);
			free(vec_data);
			free(vec_data_size);
			free(vectors);
			return(R_NilValue);				
			}
		}
	}
	
if((err=mvl_sort_indices(N, v_idx, xlength(data_list), vectors, vec_data, sort_function))!=0) {
	free(vec_data);
	free(vec_data_size);
	free(vectors);
	free(v_idx);
	error("Error sorting indices, error code %d", err);
	return(R_NilValue);		
	}
ans=PROTECT(allocVector(REALSXP, N));
pd=REAL(ans);
for(LIBMVL_OFFSET64 m=0;m<N;m++)pd[m]=v_idx[m]+1;
UNPROTECT(1);
free(vec_data);
free(vec_data_size);
free(vectors);
free(v_idx);
return(ans);	
}

SEXP hash_vectors(SEXP data_list, SEXP indices)
{
int data_idx;
LIBMVL_OFFSET64 data_offset;

double *pd;
LIBMVL_OFFSET64 *po;
int err;

void **vec_data;
LIBMVL_VECTOR **vectors;
LIBMVL_OFFSET64 *v_idx, *vec_data_length;
LIBMVL_OFFSET64 N;

SEXP ans;
	
if(TYPEOF(data_list)!=VECSXP) {
	error("order_vectors first argument must be a list of data to sort");
	return(R_NilValue);
	}

if(xlength(data_list)<1) {
	return(indices);
	}
	
if(TYPEOF(indices)!=NILSXP && xlength(indices)<1) {
	return(indices);
	}
	
vec_data=calloc(xlength(data_list), sizeof(*vec_data));
vec_data_length=calloc(xlength(data_list), sizeof(*vec_data_length));
vectors=calloc(xlength(data_list), sizeof(*vectors));
if(vec_data==NULL || vectors==NULL || vec_data_length==NULL) {
	error("Not enough memory");
	free(vec_data);
	free(vec_data_length);
	free(vectors);
	return(R_NilValue);
	}

for(LIBMVL_OFFSET64 k=0;k<xlength(data_list);k++) {
	decode_mvl_object(PROTECT(VECTOR_ELT(data_list, k)), &data_idx, &data_offset);
	UNPROTECT(1);
	vectors[k]=get_mvl_vector(data_idx, data_offset);
	
	if(vectors[k]==NULL) {
		error("Invalid MVL object in data list");
		free(vec_data);
	free(vec_data_length);
		free(vectors);
		return(R_NilValue);
		}
	vec_data[k]=libraries[data_idx].data;
	vec_data_length[k]=libraries[data_idx].length;
	}
	
if(get_indices(indices, vectors[0], &N, &v_idx)) {
	free(vec_data);
	free(vec_data_length);
	free(vectors);
	return(R_NilValue);		
	}
	
ans=PROTECT(allocVector(REALSXP, N));
pd=REAL(ans);
if((err=mvl_hash_indices(N, v_idx, (LIBMVL_OFFSET64 *)pd, xlength(data_list), vectors, vec_data, vec_data_length, LIBMVL_COMPLETE_HASH))!=0) {
	free(vec_data);
	free(vec_data_length);
	free(vectors);
	free(v_idx);
	error("Error hashing indices, code %d", err);
	UNPROTECT(1);
	return(R_NilValue);		
	}
po=(LIBMVL_OFFSET64 *)pd;
for(LIBMVL_OFFSET64 m=0;m<N;m++)po[m]=(po[m] & ((1LLU<<52)-1)) | (1023LLU<<52);
UNPROTECT(1);
free(vec_data);
free(vec_data_length);
free(vectors);
free(v_idx);
return(ans);	
}

#define N_BLOCK (1024*1024)

SEXP write_hash_vectors(SEXP idx0, SEXP data_list)
{
int data_idx, idx;
LIBMVL_OFFSET64 data_offset;

int err;

void **vec_data;
LIBMVL_VECTOR **vectors;
LIBMVL_OFFSET64 *v_idx, *hash, * vec_data_length;
LIBMVL_OFFSET64 N, offset;
double *doffset=(double *)&offset;

SEXP ans, class;
	
if(length(idx0)!=1) {
	error("fused_write_vector first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].f==NULL) {
	error("library not open for writing");
	return(R_NilValue);
	}
if(mvl_get_error(libraries[idx].ctx)!=0) {
	error("library has error status %d: %s", mvl_get_error(libraries[idx].ctx), mvl_strerror(libraries[idx].ctx));
	return(R_NilValue);
	}

if(TYPEOF(data_list)!=VECSXP) {
	error("order_vectors first argument must be a list of data to sort");
	return(R_NilValue);
	}

if(xlength(data_list)<1) {
	error("No hashes to compute");
	return(R_NilValue);
	}
	
vec_data=calloc(xlength(data_list), sizeof(*vec_data));
vec_data_length=calloc(xlength(data_list), sizeof(*vec_data_length));
vectors=calloc(xlength(data_list), sizeof(*vectors));
v_idx=calloc(N_BLOCK, sizeof(*v_idx));
hash=calloc(N_BLOCK, sizeof(*v_idx));
if(vec_data==NULL || vec_data_length==NULL || vectors==NULL || v_idx==NULL || hash==NULL) {
	error("Not enough memory");
	free(vec_data);
	free(vec_data_length);
	free(vectors);
	free(v_idx);
	free(hash);
	return(R_NilValue);
	}

for(LIBMVL_OFFSET64 k=0;k<xlength(data_list);k++) {
	decode_mvl_object(PROTECT(VECTOR_ELT(data_list, k)), &data_idx, &data_offset);
	UNPROTECT(1);
	vectors[k]=get_mvl_vector(data_idx, data_offset);
	
	if(vectors[k]==NULL) {
		error("Invalid MVL object in data list");
		free(vec_data);
		free(vec_data_length);
		free(vectors);
		free(v_idx);
		free(hash);
		return(R_NilValue);
		}
		
	vec_data[k]=libraries[data_idx].data;
	vec_data_length[k]=libraries[data_idx].length;
	}

N=mvl_vector_length(vectors[0]);

offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, N, 0, NULL, LIBMVL_NO_METADATA);
	
for(LIBMVL_OFFSET64 k=0;k<N;k+=N_BLOCK) {
	int count=N_BLOCK;
	if(k+count>mvl_vector_length(vectors[0]))count=N-k;
	for(int i=0;i<count;i++) {
		v_idx[i]=k+i;
		}
	if((err=mvl_hash_indices(count, v_idx, hash, xlength(data_list), vectors, vec_data, vec_data_length, LIBMVL_COMPLETE_HASH))!=0) {
		free(vec_data);
		free(vec_data_length);
		free(vectors);
		free(v_idx);
		free(hash);
		error("Error hashing indices, code %d", err);
		return(R_NilValue);		
		}
	mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, offset, k, count, hash);
	}
	
free(vec_data);
free(vec_data_length);
free(vectors);
free(v_idx);
free(hash);

ans=PROTECT(allocVector(REALSXP, 1));
REAL(ans)[0]=*doffset;

class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(2);

return(ans);	
}


SEXP write_groups(SEXP idx0, SEXP data_list)
{
int data_idx, idx, first_count;
LIBMVL_OFFSET64 data_offset;

int err;

void **vec_data;
LIBMVL_VECTOR **vectors;
LIBMVL_OFFSET64 *v_idx, *hash, *count, *vec_data_length;
LIBMVL_OFFSET64 N, offset;
double *doffset=(double *)&offset;
long long *first, *prev;

LIBMVL_NAMED_LIST *L;

SEXP ans, class;
	
if(length(idx0)!=1) {
	error("write_groups first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].f==NULL) {
	error("library not open for writing");
	return(R_NilValue);
	}
if(mvl_get_error(libraries[idx].ctx)!=0) {
	error("library has error status %d: %s", mvl_get_error(libraries[idx].ctx), mvl_strerror(libraries[idx].ctx));
	return(R_NilValue);
	}

if(TYPEOF(data_list)!=VECSXP) {
	error("write_groups first argument must be a list of data to sort");
	return(R_NilValue);
	}

if(xlength(data_list)<1) {
	error("No hashes to compute");
	return(R_NilValue);
	}
	
vec_data=calloc(xlength(data_list), sizeof(*vec_data));
vec_data_length=calloc(xlength(data_list), sizeof(*vec_data_length));
vectors=calloc(xlength(data_list), sizeof(*vectors));
v_idx=calloc(N_BLOCK, sizeof(*v_idx));
hash=calloc(N_BLOCK, sizeof(*v_idx));
count=calloc(N_BLOCK, sizeof(*count));
first=calloc(N_BLOCK, sizeof(*first));
prev=calloc(N_BLOCK, sizeof(*prev));
if(vec_data==NULL || vec_data_length==NULL || vectors==NULL || v_idx==NULL || hash==NULL || first==NULL || prev==NULL || count==NULL) {
	error("Not enough memory");
	free(vec_data);
	free(vec_data_length);
	free(vectors);
	free(v_idx);
	free(hash);
	free(first);
	free(prev);
	free(count);
	return(R_NilValue);
	}

for(LIBMVL_OFFSET64 k=0;k<xlength(data_list);k++) {
	decode_mvl_object(PROTECT(VECTOR_ELT(data_list, k)), &data_idx, &data_offset);
	UNPROTECT(1);
	vectors[k]=get_mvl_vector(data_idx, data_offset);
	
	if(vectors[k]==NULL) {
		error("Invalid MVL object in data list");
		free(vec_data);
		free(vec_data_length);
		free(vectors);
		free(v_idx);
		free(hash);
		free(first);
		free(prev);
		free(count);
		return(R_NilValue);
		}
		
	vec_data[k]=libraries[data_idx].data;
	vec_data_length[k]=libraries[data_idx].length;
	}

N=mvl_vector_length(vectors[0]);
if(mvl_vector_type(vectors[0])==LIBMVL_PACKED_LIST64)N--;

offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, N, 0, NULL, LIBMVL_NO_METADATA);

for(int i=0;i<N_BLOCK;i++) {
	first[i]=-1;
	count[i]=0;
	}
	
for(LIBMVL_OFFSET64 k=0;k<N;k+=N_BLOCK) {
	int bcount=N_BLOCK;
	if(k+bcount>N)bcount=N-k;
	for(int i=0;i<bcount;i++) {
		v_idx[i]=k+i;
		}
	if((err=mvl_hash_indices(bcount, v_idx, hash, xlength(data_list), vectors, vec_data, vec_data_length, LIBMVL_COMPLETE_HASH))!=0) {
		free(vec_data);
		free(vec_data_length);
		free(vectors);
		free(v_idx);
		free(hash);
		free(first);
		free(prev);
		free(count);
		error("Error hashing indices, code %d", err);
		return(R_NilValue);		
		}
	for(unsigned i=0;i<bcount;i++) {
		unsigned j=hash[i] & (N_BLOCK-1);
		count[j]++;
		if(first[j]<0) {
			first[j]=k+i;
			prev[i]=-1;
			} else {
			prev[i]=first[j]+1;
			first[j]=k+i;
			}
		}
		
	mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, offset, k, bcount, prev);
	}
	

first_count=0;
for(int i=0;i<N_BLOCK;i++) {
	if(first[i]<0)continue;
	prev[first_count]=first[i]+1;
	v_idx[first_count]=i;
	if(first_count<i)
		count[first_count]=count[i];
	first_count++;
	}
	
L=mvl_create_named_list(2);
mvl_add_list_entry(L, -1, "first", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, first_count, prev, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "mark", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, first_count, v_idx, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "count", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, first_count, count, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "prev", offset);

offset=mvl_write_named_list(libraries[idx].ctx, L);
mvl_free_named_list(L);

free(vec_data);
free(vec_data_length);
free(vectors);
free(v_idx);
free(hash);
free(first);
free(prev);
free(count);

ans=PROTECT(allocVector(REALSXP, 1));
REAL(ans)[0]=*doffset;

class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(2);

return(ans);	
}

SEXP write_spatial_groups1(SEXP idx0, SEXP data_list, SEXP sbits)
{
int data_idx, idx;
LIBMVL_OFFSET64 data_offset;
SEXP data;

void **vec_data;
LIBMVL_VECTOR **vectors;
LIBMVL_OFFSET64  *count;
double *values;
LIBMVL_OFFSET64 N, offset;
double *doffset=(double *)&offset;
long long *first, *prev;
int *hash, *bits;
LIBMVL_VEC_STATS *vstats;
unsigned block_length, block_bits;

LIBMVL_NAMED_LIST *L;

SEXP ans, class;
	
if(length(idx0)!=1) {
	error("fused_write_vector first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].f==NULL) {
	error("library not open for writing");
	return(R_NilValue);
	}
if(mvl_get_error(libraries[idx].ctx)!=0) {
	error("library has error status %d: %s", mvl_get_error(libraries[idx].ctx), mvl_strerror(libraries[idx].ctx));
	return(R_NilValue);
	}

if(TYPEOF(data_list)!=VECSXP) {
	error("order_vectors first argument must be a list of data to sort");
	return(R_NilValue);
	}

if(xlength(data_list)<1) {
	error("No hashes to compute");
	return(R_NilValue);
	}
	
if(xlength(data_list)!=xlength(sbits)) {
	error("Need to specify number of useful bits for each vector");
	return(R_NilValue);
	}
	
bits=INTEGER(sbits);
block_bits=0;
for(LIBMVL_OFFSET64 i=0; i< xlength(data_list);i++)block_bits+=bits[i];
if(block_bits>30) {
	error("Too many bits: %d total", block_bits);
	return(R_NilValue);
	}
block_length=1<<block_bits;
	
vec_data=calloc(xlength(data_list), sizeof(*vec_data));
vectors=calloc(xlength(data_list), sizeof(*vectors));
vstats=calloc(xlength(data_list), sizeof(*vstats));
values=calloc(N_BLOCK, sizeof(*values));
hash=calloc(N_BLOCK, sizeof(*hash));
count=calloc(block_length, sizeof(*count));
first=calloc(block_length, sizeof(*first));
prev=calloc(N_BLOCK, sizeof(*prev));
if(vec_data==NULL || vectors==NULL || vstats==NULL || hash==NULL || first==NULL || prev==NULL || values==NULL) {
	error("Not enough memory");
	free(vec_data);
	free(vectors);
	free(vstats);
	free(hash);
	free(first);
	free(prev);
	free(values);
	return(R_NilValue);
	}

for(LIBMVL_OFFSET64 k=0;k<xlength(data_list);k++) {
	decode_mvl_object(PROTECT(VECTOR_ELT(data_list, k)), &data_idx, &data_offset);
	UNPROTECT(1);
	vectors[k]=get_mvl_vector(data_idx, data_offset);
	
	if(vectors[k]==NULL) {
		error("Invalid MVL object in data list");
		free(vec_data);
		free(vectors);
		free(hash);
		free(first);
		free(prev);
		free(count);
		free(values);
		return(R_NilValue);
		}
	vec_data[k]=libraries[data_idx].data;
	mvl_compute_vec_stats(vectors[k], &(vstats[k]));
	}

N=mvl_vector_length(vectors[0]);
for(LIBMVL_OFFSET64 k=1;k<xlength(data_list);k++) {
	if(mvl_vector_length(vectors[k])!=N) {
		error("All MVL vectors should be equal length");
		free(vec_data);
		free(vectors);
		free(hash);
		free(first);
		free(prev);
		free(count);
		free(values);
		return(R_NilValue);
		}
	}

offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, N, 0, NULL, LIBMVL_NO_METADATA);

for(int i=0;i<block_length;i++) {
	first[i]=-1;
	count[i]=0;
	}
	
for(LIBMVL_OFFSET64 k=0;k<N;k+=N_BLOCK) {
	int bcount=N_BLOCK;
	if(k+bcount>N)bcount=N-k;
	
	memset(hash, 0, sizeof(*hash)*bcount);
	
	for(LIBMVL_OFFSET64 m=0;m<xlength(data_list);m++) {
		mvl_normalize_vector(vectors[m], &(vstats[m]), k, k+bcount, values);
		int shift=bits[m];
		int mult=1<<shift;
		int mask=mult-1;
		for(int i=0;i<bcount;i++) {
			hash[i]=hash[i]<<shift | ((int)floor(values[i]*mult) & mask);
			}
		}
	
	for(unsigned i=0;i<bcount;i++) {
		unsigned j=hash[i] & (block_length-1);
		count[j]++;
		if(first[j]<0) {
			first[j]=k+i;
			prev[i]=-1;
			} else {
			prev[i]=first[j]+1;
			first[j]=k+i;
			}
		}
		
	mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, offset, k, bcount, prev);
	}
	

// first_count=0;
// for(int i=0;i<block_length;i++) {
// 	if(first[i]<0)continue;
// 	first[first_count]=first[i]+1;
// 	v_idx[first_count]=i;
// 	if(first_count<i)
// 		count[first_count]=count[i];
// 	first_count++;
// 	}
	
L=mvl_create_named_list(2);
mvl_add_list_entry(L, -1, "first", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, block_length, first, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "count", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, block_length, count, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "bits", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT32, xlength(data_list), bits, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "vector_stats", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_DOUBLE, xlength(data_list)*sizeof(*vstats)/sizeof(double), vstats, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "prev", offset);

offset=mvl_write_named_list(libraries[idx].ctx, L);
mvl_free_named_list(L);

free(vec_data);
free(vectors);
free(hash);
free(first);
free(prev);
free(count);
free(values);

ans=PROTECT(allocVector(REALSXP, 1));
REAL(ans)[0]=*doffset;

class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(2);

return(ans);	
}


SEXP write_spatial_groups(SEXP idx0, SEXP data_list, SEXP sbits)
{
int data_idx, idx;
LIBMVL_OFFSET64 data_offset;

void **vec_data;
LIBMVL_VECTOR **vectors;
LIBMVL_OFFSET64  *count, *hash, *first_hash, max_count;
double *values;
LIBMVL_OFFSET64 N, offset, first_idx_size, first_idx_free, N2;
double *doffset=(double *)&offset;
long long *first_idx, *prev, *first_idx_prev, *first;
int *bits;
LIBMVL_VEC_STATS *vstats;
unsigned block_length, block_bits;

LIBMVL_NAMED_LIST *L;

SEXP ans, class;
	
if(length(idx0)!=1) {
	error("fused_write_vector first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].f==NULL) {
	error("library not open for writing");
	return(R_NilValue);
	}
if(mvl_get_error(libraries[idx].ctx)!=0) {
	error("library has error status %d: %s", mvl_get_error(libraries[idx].ctx), mvl_strerror(libraries[idx].ctx));
	return(R_NilValue);
	}

if(TYPEOF(data_list)!=VECSXP) {
	error("order_vectors first argument must be a list of data to sort");
	return(R_NilValue);
	}

if(xlength(data_list)<1) {
	error("No hashes to compute");
	return(R_NilValue);
	}
	
if(xlength(data_list)!=xlength(sbits)) {
	error("Need to specify number of useful bits for each vector");
	return(R_NilValue);
	}
	
bits=INTEGER(sbits);
block_bits=0;
for(LIBMVL_OFFSET64 i=0; i< xlength(data_list);i++)block_bits+=bits[i];
if(block_bits>63) {
	error("Too many bits: %d total", block_bits);
	return(R_NilValue);
	}
	
if(block_bits>25)block_length=1<<25;
	else
	block_length=1<<block_bits;

first_idx_free=0;
first_idx_size=N_BLOCK;

first=calloc(first_idx_size, sizeof(*first_idx));
first_idx_prev=calloc(first_idx_size, sizeof(*first_idx_prev));
first_hash=calloc(first_idx_size, sizeof(*first_hash));
count=calloc(first_idx_size, sizeof(*count));
	
vec_data=calloc(xlength(data_list), sizeof(*vec_data));
vectors=calloc(xlength(data_list), sizeof(*vectors));
vstats=calloc(xlength(data_list), sizeof(*vstats));
values=calloc(N_BLOCK, sizeof(*values));
hash=calloc(N_BLOCK, sizeof(*hash));
first_idx=calloc(block_length, sizeof(*first));
prev=calloc(N_BLOCK, sizeof(*prev));
if(vec_data==NULL || vectors==NULL || vstats==NULL || hash==NULL || first==NULL || prev==NULL || values==NULL || first_idx==NULL || first_idx_prev==NULL || first_hash==NULL) {
	error("Not enough memory");
	free(vec_data);
	free(vectors);
	free(vstats);
	free(hash);
	free(first);
	free(prev);
	free(values);
	free(first_idx);
	free(first_idx_prev);
	free(first_hash);
	return(R_NilValue);
	}

for(LIBMVL_OFFSET64 k=0;k<xlength(data_list);k++) {
	decode_mvl_object(PROTECT(VECTOR_ELT(data_list, k)), &data_idx, &data_offset);
	UNPROTECT(1);
	vectors[k]=get_mvl_vector(data_idx, data_offset);
	
	if(vectors[k]==NULL) {
		error("Invalid MVL object in data list");
		free(vec_data);
		free(vectors);
		free(hash);
		free(first);
		free(prev);
		free(count);
		free(values);
		free(first_idx);
		free(first_idx_prev);
		free(first_hash);
		return(R_NilValue);
		}
	vec_data[k]=libraries[data_idx].data;
	mvl_compute_vec_stats(vectors[k], &(vstats[k]));
	}

N=mvl_vector_length(vectors[0]);
for(LIBMVL_OFFSET64 k=1;k<xlength(data_list);k++) {
	if(mvl_vector_length(vectors[k])!=N) {
		error("All MVL vectors should be equal length");
		free(vec_data);
		free(vectors);
		free(hash);
		free(first);
		free(prev);
		free(count);
		free(values);
		free(first_idx);
		free(first_idx_prev);
		free(first_hash);
		return(R_NilValue);
		}
	}

offset=mvl_start_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, N, 0, NULL, LIBMVL_NO_METADATA);

for(int i=0;i<block_length;i++) {
	first_idx[i]=-1;
	}
	
for(LIBMVL_OFFSET64 k=0;k<N;k+=N_BLOCK) {
	int bcount=N_BLOCK;
	if(k+bcount>N)bcount=N-k;
	
	memset(hash, 0, sizeof(*hash)*bcount);
	
	for(LIBMVL_OFFSET64 m=0;m<xlength(data_list);m++) {
		mvl_normalize_vector(vectors[m], &(vstats[m]), k, k+bcount, values);
		int shift=bits[m];
		int mult=1<<shift;
		int mask=mult-1;
		for(int i=0;i<bcount;i++) {
			int m0=floor(values[i]*mult)-mult;
			if(m0<0)m0=0;
			if(m0>mask)m0=mask;
			hash[i]=hash[i]<<shift | (m0);
			}
		}
	
	for(unsigned i=0;i<bcount;i++) {
		unsigned j=mvl_randomize_bits64(hash[i]) & (block_length-1);
		long long idx;
		idx=first_idx[j];
		while(idx>=0) {
			if(first_hash[idx]==hash[i])break;
			idx=first_idx_prev[idx];
			}
		if(idx<0) {
			idx=first_idx_free;
			if(first_idx_free>=first_idx_size) {
				first_idx_size=2*first_idx_size+256;
				LIBMVL_OFFSET64 *p;
				
				p=calloc(first_idx_size, sizeof(*p));
				if(p==NULL)goto cleanup1;
				memcpy(p, first_hash, first_idx_free*sizeof(*p));
				free(first_hash);
				first_hash=p;
				
				p=calloc(first_idx_size, sizeof(*p));
				if(p==NULL)goto cleanup1;
				memcpy(p, count, first_idx_free*sizeof(*p));
				free(count);
				count=p;

				p=calloc(first_idx_size, sizeof(*p));
				if(p==NULL)goto cleanup1;
				memcpy(p, first_idx_prev, first_idx_free*sizeof(*p));
				free(first_idx_prev);
				first_idx_prev=(long long *)p;
				
				p=calloc(first_idx_size, sizeof(*p));
				if(p==NULL)goto cleanup1;
				memcpy(p, first, first_idx_free*sizeof(*p));
				free(first);
				first=(long long *)p;
				}
			first_idx_free++;
			count[idx]=0;
			first_hash[idx]=hash[i];
			first_idx_prev[idx]=first_idx[j];
			first_idx[j]=idx;
			first[idx]=-1;
			}
			
		count[idx]++;
		
		
		if(first[idx]<0) {
			first[idx]=k+i+1;
			prev[i]=-1;
			} else {
			prev[i]=first[idx];
			first[idx]=k+i+1;
			}
		}
		
	mvl_rewrite_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, offset, k, bcount, prev);
	}
	

// first_count=0;
// for(int i=0;i<block_length;i++) {
// 	if(first[i]<0)continue;
// 	first[first_count]=first[i]+1;
// 	v_idx[first_count]=i;
// 	if(first_count<i)
// 		count[first_count]=count[i];
// 	first_count++;
// 	}
	
L=mvl_create_named_list(2);
mvl_add_list_entry(L, -1, "index_type", MVL_WVEC(libraries[idx].ctx, LIBMVL_VECTOR_INT32, MVL_SPATIAL_INDEX1));
mvl_add_list_entry(L, -1, "first", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, first_idx_free, first, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "count", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, first_idx_free, count, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "mark", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, first_idx_free, first_hash, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "bits", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT32, xlength(data_list), bits, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "vector_stats", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_DOUBLE, xlength(data_list)*sizeof(*vstats)/sizeof(double), vstats, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "prev", offset);

max_count=0;
for(LIBMVL_OFFSET64 i=0;i<first_idx_free;i++)
	if(count[i]>max_count)max_count=count[i];

mvl_add_list_entry(L, -1, "max_count", MVL_WVEC(libraries[idx].ctx, LIBMVL_VECTOR_INT64, max_count));
	
N2=1;
while(N2<first_idx_free)N2=N2<<1;
if(N2>first_idx_size) {
	error("Internal error: N2>first_idx_size (%lld > %lld)\n", N2, first_idx_size);
	} else {
	for(LIBMVL_OFFSET64 i=0;i<N2;i++)first[i]=-1;
	
	for(LIBMVL_OFFSET64 i=0;i<first_idx_free;i++) {
		LIBMVL_OFFSET64 j=mvl_randomize_bits64(first_hash[i]) & (N2-1);
		
		if(first[j]<0) {
			first_idx_prev[i]=-1;
			first[j]=i+1;
			} else {
			first_idx_prev[i]=first[j];
			first[j]=i+1;
			}
		}
	}
mvl_add_list_entry(L, -1, "first_mark", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, N2, first, LIBMVL_NO_METADATA));
mvl_add_list_entry(L, -1, "prev_mark", mvl_write_vector(libraries[idx].ctx, LIBMVL_VECTOR_INT64, first_idx_free, first_idx_prev, LIBMVL_NO_METADATA));

offset=mvl_write_named_list2(libraries[idx].ctx, L, "MVL_INDEX");
mvl_free_named_list(L);

free(vec_data);
free(vectors);
free(hash);
free(first);
free(prev);
free(count);
free(values);
free(first_idx);
free(first_idx_prev);
free(first_hash);

ans=PROTECT(allocVector(REALSXP, 1));
REAL(ans)[0]=*doffset;

class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(2);

return(ans);	

cleanup1: {
	error("Not enough memory");
	free(vec_data);
	free(vectors);
	free(vstats);
	free(hash);
	free(first);
	free(prev);
	free(values);
	free(first_idx);
	free(first_idx_prev);
	free(first_hash);
	return(R_NilValue);
	}
}

void normalize_vector(SEXP sexp, LIBMVL_VEC_STATS *stats, LIBMVL_OFFSET64 i0, LIBMVL_OFFSET64 i1, double *out)
{
int data_idx;
LIBMVL_OFFSET64 data_offset;
LIBMVL_VECTOR *vec;
double *pd;
int *pi;
double scale, center;
scale=0.5*stats->scale;
center=1.5-stats->center*scale;
//fprintf(stderr, "scale=%g center=%g\n", scale, center);

if(i0>=i1)return;
switch(TYPEOF(sexp)) {
	case VECSXP: 
		decode_mvl_object(sexp, &data_idx, &data_offset);
		vec=get_mvl_vector(data_idx, data_offset);
		if(vec==NULL) {
			error("Provided vector is a list and not an MVL object");
			memset(out, 0, (i1-i0)*sizeof(*out));
			return;
			}
		mvl_normalize_vector(vec, stats, i0, i1, out);
		return;
	case REALSXP:
		pd=REAL(sexp);
		if(i1>xlength(sexp)) {
			error("Vector lengths do not match");
			memset(out, 0, (i1-i0)*sizeof(*out));
			i1=xlength(sexp);
			}
		for(LIBMVL_OFFSET64 i=i0;i<i1;i++)
			out[i-i0]=pd[i]*scale+center;
		return;
	case INTSXP:
		pi=INTEGER(sexp);
		if(i1>xlength(sexp)) {
			error("Vector lengths do not match");
			memset(out, 0, (i1-i0)*sizeof(*out));
			i1=xlength(sexp);
			}
		for(LIBMVL_OFFSET64 i=i0;i<i1;i++)
			out[i-i0]=pi[i]*scale+center;
		return;
	default:
		error("Cannot handle R vector of type %d", TYPEOF(sexp));
		memset(out, 0, (i1-i0)*sizeof(*out));
		return;
	}
}

/* Pass -1 to expected type to disable it */
int hash_vector_range(SEXP sexp, LIBMVL_OFFSET64 i0, LIBMVL_OFFSET64 i1, int expected_type, LIBMVL_OFFSET64 *out)
{
int data_idx;
LIBMVL_OFFSET64 data_offset;
LIBMVL_VECTOR *vec;
double *pd;
int *pi;
unsigned char *ps;
int err;
long long da;
double dd;
int warn_once=1;
SEXP sch;

if(i0>=i1)return 0;

switch(TYPEOF(sexp)) {
	case VECSXP: 
		decode_mvl_object(sexp, &data_idx, &data_offset);
		vec=get_mvl_vector(data_idx, data_offset);
		if(vec==NULL) {
			error("Provided vector is a list and not an MVL object");
			memset(out, 0, (i1-i0)*sizeof(*out));
			return -1;
			}
		if(0 && mvl_vector_type(vec)!=expected_type) {
			error("Vector types do not match: expected %d got %d", mvl_vector_type(vec), expected_type);
			}
		if((err=mvl_hash_range(i0, i1, out, 1, &vec, (void **)&(libraries[data_idx].data), &(libraries[data_idx].length), LIBMVL_ACCUMULATE_HASH))!=0) {
			error("Error computing hashes (%d)", err);
			memset(out, 0, (i1-i0)*sizeof(*out));
			return -1;
			}
		return 0;
	case REALSXP:
		pd=REAL(sexp);
		if(i1>xlength(sexp)) {
			error("Vector lengths do not match");
			memset(out, 0, (i1-i0)*sizeof(*out));
			i1=xlength(sexp);
			}
		switch(expected_type) {
			case LIBMVL_VECTOR_INT32:
			case LIBMVL_VECTOR_INT64:
				for(LIBMVL_OFFSET64 i=i0;i<i1;i++) {
					dd=floor(pd[i]);
					if(dd!=pd[i] && warn_once) {
						error("numeric() values are not integer when quering integer vector");
						warn_once=0;
						}
					da=dd;
					out[i-i0]=mvl_accumulate_int64_hash64(out[i-i0], &(da), 1);
					}
				break;
			case LIBMVL_VECTOR_FLOAT:
			case LIBMVL_VECTOR_DOUBLE:
				for(LIBMVL_OFFSET64 i=i0;i<i1;i++)
					out[i-i0]=mvl_accumulate_double_hash64(out[i-i0], &(pd[i]), 1);
				break;
			default:
				if(expected_type>=0) {
					error("using numeric() values to query MVL vector of type %d", expected_type);
					warn_once=0;
					}
				for(LIBMVL_OFFSET64 i=i0;i<i1;i++)
					out[i-i0]=mvl_accumulate_double_hash64(out[i-i0], &(pd[i]), 1);
				break;
			}
		return 0;
	case INTSXP:
		pi=INTEGER(sexp);
		if(i1>xlength(sexp)) {
			error("Vector lengths do not match");
			memset(out, 0, (i1-i0)*sizeof(*out));
			i1=xlength(sexp);
			}
		switch(expected_type) {
			case LIBMVL_VECTOR_INT32:
			case LIBMVL_VECTOR_INT64:
				for(LIBMVL_OFFSET64 i=i0;i<i1;i++)
					out[i-i0]=mvl_accumulate_int32_hash64(out[i-i0], &(pi[i]), 1);
				break;
			case LIBMVL_VECTOR_FLOAT:
			case LIBMVL_VECTOR_DOUBLE:
				for(LIBMVL_OFFSET64 i=i0;i<i1;i++) {
					dd=pi[i];
					out[i-i0]=mvl_accumulate_double_hash64(out[i-i0], &(dd), 1);
					}
				break;
			default:
				if(expected_type>=0) {
					error("using numeric() values to query MVL vector of type %d", expected_type);
					warn_once=0;
					}
				for(LIBMVL_OFFSET64 i=i0;i<i1;i++)
					out[i-i0]=mvl_accumulate_int32_hash64(out[i-i0], &(pi[i]), 1);
				break;
			}
		return 0;
	case STRSXP:
		if(i1>xlength(sexp)) {
			error("Vector lengths do not match");
			memset(out, 0, (i1-i0)*sizeof(*out));
			i1=xlength(sexp);
			}
		/* TODO: NA_STRING in R is just "NA" we might, or might not want to alter the code below */
		for(LIBMVL_OFFSET64 i=i0;i<i1;i++) {
			sch=STRING_ELT(sexp, i);
			if(sch==NA_STRING) 
				out[i-i0]=mvl_accumulate_hash64(out[i-i0], (unsigned char*)MVL_NA_STRING, MVL_NA_STRING_LENGTH);
				else {
				ps=(unsigned char *)CHAR(sch);
				out[i-i0]=mvl_accumulate_hash64(out[i-i0], ps, strlen((char *)ps));
				}
			}
		
		return 0;
	default:
		error("Cannot handle R vector of type %d", TYPEOF(sexp));
		memset(out, 0, (i1-i0)*sizeof(*out));
		return 0;
	}
}

SEXP get_neighbors(SEXP spatial_index, SEXP data_list) 
{
LIBMVL_OFFSET64 data_offset, index_offset, *query_mark, Nv, N2, indices_size, indices_free, *indices, ball_size, Nbits;
int data_idx, index_idx;
LIBMVL_VECTOR *vec_bits, *vec_first, *vec_first_mark, *vec_prev_mark, *vec_mark, *vec_prev, *vec, *vec_stats, *vec_max_count;
SEXP ans, sa;
double *values, *pd;
LIBMVL_NAMED_LIST *L;
int *bits;
long long *first, *first_mark, *prev_mark, *mark, *prev, max_count;
LIBMVL_VEC_STATS *vstats;
char *ball;

if(TYPEOF(data_list)!=VECSXP) {
	error("Second argument should be a list (or data frame) of vectors to query");
	return(R_NilValue);
	}

decode_mvl_object(spatial_index, &index_idx, &index_offset);
L=get_mvl_named_list(index_idx, index_offset);

if(L==NULL) {
	error("Not a spatial index (1)");
	return(R_NilValue);
	}

vec_bits=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "bits"));
vec_mark=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "mark"));
vec_first=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "first"));
vec_first_mark=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "first_mark"));
vec_prev_mark=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "prev_mark"));
vec_prev=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "prev"));
vec_stats=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "vector_stats"));
vec_max_count=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "max_count"));

if(vec_bits==NULL || vec_mark==NULL || vec_first==NULL || vec_first_mark==NULL || vec_prev==NULL || vec_stats==NULL || vec_max_count==NULL) {
	mvl_free_named_list(L);
	error("Not a spatial index (2)");
	return(R_NilValue);
	}

if(mvl_vector_type(vec_bits)!=LIBMVL_VECTOR_INT32 
	|| mvl_vector_type(vec_first)!=LIBMVL_VECTOR_INT64
	|| mvl_vector_type(vec_mark)!=LIBMVL_VECTOR_INT64
	|| mvl_vector_type(vec_first_mark)!=LIBMVL_VECTOR_INT64
	|| mvl_vector_type(vec_prev_mark)!=LIBMVL_VECTOR_INT64
	|| mvl_vector_type(vec_prev)!=LIBMVL_VECTOR_INT64
	|| mvl_vector_type(vec_stats)!=LIBMVL_VECTOR_DOUBLE
	|| mvl_vector_type(vec_max_count)!=LIBMVL_VECTOR_INT64
	) {
	mvl_free_named_list(L);
	error("Not a spatial index (3)");
	return(R_NilValue);
	}
	
N2=mvl_vector_length(vec_first_mark);
	
bits=mvl_vector_data_int32(vec_bits);
Nbits=mvl_vector_length(vec_bits);

if(Nbits!=length(data_list) || length(data_list)==0) {
	mvl_free_named_list(L);
	error("Query columns do not match spatial index: %d vs %llu", length(data_list), Nbits);
	return(R_NilValue);
	}
	
if(Nbits*sizeof(LIBMVL_VEC_STATS)!=mvl_vector_length(vec_stats)*sizeof(double)) {
	mvl_free_named_list(L);
	error("Not a spatial index (4) : malformed vector statistics");
	return(R_NilValue);
	}
	
/* This cutoff is determined by the call to calloc() to allocate ball below */
if(Nbits>36) {
	mvl_free_named_list(L);
	error("Nbits=%llu is too large (>36)", Nbits);
	return(R_NilValue);
	}
	
	
vstats=(LIBMVL_VEC_STATS *)mvl_vector_data_double(vec_stats);

first_mark=mvl_vector_data_int64(vec_first_mark);
prev_mark=mvl_vector_data_int64(vec_prev_mark);
mark=mvl_vector_data_int64(vec_mark);
prev=mvl_vector_data_int64(vec_prev);
first=mvl_vector_data_int64(vec_first);
max_count=mvl_vector_data_int64(vec_max_count)[0];


switch(TYPEOF(VECTOR_ELT(data_list, 0))) {
	case REALSXP:
	case INTSXP:
		Nv=xlength(VECTOR_ELT(data_list, 0));
		break;
	case VECSXP:
		decode_mvl_object(PROTECT(VECTOR_ELT(data_list, 0)), &data_idx, &data_offset);
		UNPROTECT(1);
		vec=get_mvl_vector(data_idx, data_offset);
		if(vec==NULL) {
			mvl_free_named_list(L);
			error("Not an MVL object");
			return(R_NilValue);
			}
		Nv=mvl_vector_length(vec);
		if(mvl_vector_type(vec)==LIBMVL_PACKED_LIST64)Nv--;
		break;
	default:
		mvl_free_named_list(L);
		error("Cannot handle R vector of type %d", TYPEOF(VECTOR_ELT(data_list, 0)));
		return(R_NilValue);
		break;
	}
	
ball_size=1;
for(int i=0;i<Nbits;i++)
	ball_size*=3;

values=calloc(Nv, sizeof(*values));
query_mark=calloc(Nv, sizeof(*query_mark));
indices_size=max_count*ball_size;
indices=calloc(indices_size, sizeof(*indices));
ball=calloc(ball_size*Nbits, sizeof(*ball));

if(values==NULL || query_mark==NULL || indices==NULL || ball==NULL) {
	error("Not enough memory");
	free(values);
	free(query_mark);
	free(indices);
	free(ball);
	mvl_free_named_list(L);
	return(R_NilValue);
	}
	
for(long long b=0;b<ball_size; b++) {
	long long b0=b;
	for(int i=0;i<Nbits;i++) {
		ball[b*Nbits+i]=b0 % 3;
		b0=b0/3;
		}
	}
	
memset(query_mark, 0, Nv*sizeof(*query_mark));

for(int i=0;i<Nbits;i++) {
	int shift=bits[i];
	int mult=1<<shift;
	int mask=mult-1;
	normalize_vector(PROTECT(VECTOR_ELT(data_list, i)), &(vstats[i]), 0, Nv, values);
	UNPROTECT(1);
	for(LIBMVL_OFFSET64 k=0;k<Nv;k++) {
		int m=floor(values[k]*mult)-mult;
		if(m<0)m=0;
		if(m>mask)m=mask;
		query_mark[k]=query_mark[k]<<shift | (m);
		//fprintf(stderr, "%d %lld m=%d shift=%d mult=%d mask=%d mark=%lld value=%g\n", i, k, m, shift, mult, mask, query_mark[k], values[k]);
		}
	}

ans=PROTECT(allocVector(VECSXP, Nv));	
	
for(LIBMVL_OFFSET64 i=0;i<Nv;i++) {
	indices_free=0;
	
	for(long long b=0;b<ball_size;b++) {
		LIBMVL_OFFSET64 center_mark=query_mark[i];
		LIBMVL_OFFSET64 nmark=0;
		int shift=0;
		int skip=0;
		
		for(int j=Nbits-1;j>=0;j--) {
			LIBMVL_OFFSET64 mask=(1<<bits[j])-1;
			LIBMVL_OFFSET64 pm=(center_mark>>shift) & mask;
			char nudge=ball[b*Nbits+j];
			if(nudge==2) {
				if(pm==mask) {
					skip=1;
					break;
					}
				pm++;
				} else
			if(nudge==0) {
				if(pm==0) {
					skip=1;
					break;
					}
				pm--;
				}
			nmark|=pm<<shift;
			shift+=bits[j];
			}
		if(skip)continue;
			
		//fprintf(stderr, "b=%lld center=0x%08llx mark=%08llx (%lld %lld)\n", b, center_mark, nmark, center_mark, nmark);
		
		LIBMVL_OFFSET64 j=mvl_randomize_bits64(nmark) & (N2-1);
		long long k;

		
		//fprintf(stderr, "%lld %lld\n", query_mark[i], j);
		
		k=first_mark[j];
		while(k>=0) {
			if(mark[k-1]==nmark)break;
			k=prev_mark[k-1];
			}
		if(k<0)continue;
		
	//	fprintf(stderr, "%lld %lld %lld\n", query_mark[i], j, k);
		
		for(long long m=first[k-1]; m>=0 ; m=prev[m-1]) {
			indices[indices_free]=m;
			indices_free++;
			}
			
		if(indices_free>indices_size) {
			Rprintf("*** INTERNAL ERROR: array overflow");
			}
		}
	
	sa=PROTECT(allocVector(REALSXP, indices_free));
	pd=REAL(sa);
	for(LIBMVL_OFFSET64 m=0;m<indices_free;m++)
		pd[m]=indices[m];
	SET_VECTOR_ELT(ans, i, sa);
	UNPROTECT(1);
	}
	
free(values);
free(query_mark);
free(indices);
free(ball);
mvl_free_named_list(L);

UNPROTECT(1);
return(ans);
}

SEXP neighbors_lapply(SEXP spatial_index, SEXP data_list, SEXP fn, SEXP env) 
{
LIBMVL_OFFSET64 data_offset, index_offset, *query_mark, Nv, N2, indices_size, indices_free, *indices, Nbits, ball_size;
int data_idx, index_idx;
LIBMVL_VECTOR *vec_bits, *vec_first, *vec_first_mark, *vec_prev_mark, *vec_mark, *vec_prev, *vec, *vec_stats, *vec_max_count;
SEXP ans, sa, sa2, R_fcall, tmp;
double *values, *pd, *pd2;
LIBMVL_NAMED_LIST *L;
int *bits;
long long *first, *first_mark, *prev_mark, *mark, *prev, max_count;
LIBMVL_VEC_STATS *vstats;
char *ball;

if(TYPEOF(data_list)!=VECSXP) {
	error("Second argument should be a list (or data frame) of vectors to query");
	return(R_NilValue);
	}

decode_mvl_object(spatial_index, &index_idx, &index_offset);
L=get_mvl_named_list(index_idx, index_offset);

if(L==NULL) {
	error("Not a spatial index (1)");
	return(R_NilValue);
	}

vec_bits=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "bits"));
vec_mark=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "mark"));
vec_first=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "first"));
vec_first_mark=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "first_mark"));
vec_prev_mark=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "prev_mark"));
vec_prev=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "prev"));
vec_stats=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "vector_stats"));
vec_max_count=get_mvl_vector(index_idx, mvl_find_list_entry(L, -1, "max_count"));

if(vec_bits==NULL || vec_mark==NULL || vec_first==NULL || vec_first_mark==NULL || vec_prev==NULL || vec_stats==NULL || vec_max_count==NULL) {
	mvl_free_named_list(L);
	error("Not a spatial index (2)");
	return(R_NilValue);
	}

if(mvl_vector_type(vec_bits)!=LIBMVL_VECTOR_INT32 
	|| mvl_vector_type(vec_first)!=LIBMVL_VECTOR_INT64
	|| mvl_vector_type(vec_mark)!=LIBMVL_VECTOR_INT64
	|| mvl_vector_type(vec_first_mark)!=LIBMVL_VECTOR_INT64
	|| mvl_vector_type(vec_prev_mark)!=LIBMVL_VECTOR_INT64
	|| mvl_vector_type(vec_prev)!=LIBMVL_VECTOR_INT64
	|| mvl_vector_type(vec_stats)!=LIBMVL_VECTOR_DOUBLE
	|| mvl_vector_type(vec_max_count)!=LIBMVL_VECTOR_INT64
	) {
	mvl_free_named_list(L);
	error("Not a spatial index (3)");
	return(R_NilValue);
	}
	
N2=mvl_vector_length(vec_first_mark);
	
bits=mvl_vector_data_int32(vec_bits);
Nbits=mvl_vector_length(vec_bits);

if(Nbits!=length(data_list) || length(data_list)==0) {
	mvl_free_named_list(L);
	error("Query columns do not match spatial index: %d vs %llu", length(data_list), Nbits);
	return(R_NilValue);
	}
	
if(Nbits*sizeof(LIBMVL_VEC_STATS)!=mvl_vector_length(vec_stats)*sizeof(double)) {
	mvl_free_named_list(L);
	error("Not a spatial index (4) : malformed vectors statistics");
	return(R_NilValue);
	}
	
vstats=(LIBMVL_VEC_STATS *)mvl_vector_data_float(vec_stats);

first_mark=mvl_vector_data_int64(vec_first_mark);
prev_mark=mvl_vector_data_int64(vec_prev_mark);
mark=mvl_vector_data_int64(vec_mark);
prev=mvl_vector_data_int64(vec_prev);
first=mvl_vector_data_int64(vec_first);
max_count=mvl_vector_data_int64(vec_max_count)[0];


switch(TYPEOF(VECTOR_ELT(data_list, 0))) {
	case REALSXP:
	case INTSXP:
		Nv=xlength(VECTOR_ELT(data_list, 0));
		break;
	case VECSXP:
		decode_mvl_object(PROTECT(VECTOR_ELT(data_list, 0)), &data_idx, &data_offset);
		UNPROTECT(1);
		vec=get_mvl_vector(data_idx, data_offset);
		if(vec==NULL) {
			mvl_free_named_list(L);
			error("Not an MVL object");
			return(R_NilValue);
			}
		Nv=mvl_vector_length(vec);
		if(mvl_vector_type(vec)==LIBMVL_PACKED_LIST64)Nv--;
		break;
	default:
		mvl_free_named_list(L);
		error("Cannot handle R vector of type %d", TYPEOF(VECTOR_ELT(data_list, 0)));
		return(R_NilValue);
		break;
	}
	
ball_size=1;
for(int i=0;i<Nbits;i++)
	ball_size*=3;

values=calloc(Nv, sizeof(*values));
query_mark=calloc(Nv, sizeof(*query_mark));
indices_size=max_count*ball_size;
indices=calloc(indices_size, sizeof(*indices));
ball=calloc(ball_size*Nbits, sizeof(*ball));
for(long long b=0;b<ball_size; b++) {
	long long b0=b;
	for(int i=0;i<Nbits;i++) {
		ball[b*Nbits+i]=b0 % 3;
		b0=b0/3;
		}
	}

if(values==NULL || query_mark==NULL || indices==NULL || ball==NULL) {
	error("Not enough memory");
	free(values);
	free(query_mark);
	free(indices);
	free(ball);
	mvl_free_named_list(L);
	return(R_NilValue);
	}
	
memset(query_mark, 0, Nv*sizeof(*query_mark));

for(int i=0;i<Nbits;i++) {
	int shift=bits[i];
	int mult=1<<shift;
	int mask=mult-1;
	normalize_vector(PROTECT(VECTOR_ELT(data_list, i)), &(vstats[i]), 0, Nv, values);
	UNPROTECT(1);
	for(LIBMVL_OFFSET64 k=0;k<Nv;k++) {
		int m=floor(values[k]*mult)-mult;
		if(m<0)m=0;
		if(m>mask)m=mask;
		query_mark[k]=query_mark[k]<<shift | (m);
		//fprintf(stderr, "%d %lld m=%d shift=%d mult=%d mask=%d mark=%lld value=%g\n", i, k, m, shift, mult, mask, query_mark[k], values[k]);
		}
	}

ans=PROTECT(allocVector(VECSXP, Nv));	

R_fcall = PROTECT(lang3(fn, R_NilValue, R_NilValue));
	
sa=PROTECT(allocVector(REALSXP, indices_size));
// ENABLE_REFCNT(sa);
// INCREMENT_REFCNT(sa);
// INCREMENT_NAMED(sa);
pd=REAL(sa);



for(LIBMVL_OFFSET64 i=0;i<Nv;i++) {
	indices_free=0;
	
	for(long long b=0;b<ball_size;b++) {
		LIBMVL_OFFSET64 center_mark=query_mark[i];
		LIBMVL_OFFSET64 nmark=0;
		int shift=0;
		int skip=0;
		
		for(int j=Nbits-1;j>=0;j--) {
			LIBMVL_OFFSET64 mask=(1<<bits[j])-1;
			LIBMVL_OFFSET64 pm=(center_mark>>shift) & mask;
			char nudge=ball[b*Nbits+j];
			if(nudge==2) {
				if(pm==mask) {
					skip=1;
					break;
					}
				pm++;
				} else
			if(nudge==0) {
				if(pm==0) {
					skip=1;
					break;
					}
				pm--;
				}
			nmark|=pm<<shift;
			shift+=bits[j];
			}
		if(skip)continue;
			
		//fprintf(stderr, "b=%lld center=0x%08llx mark=%08llx (%lld %lld)\n", b, center_mark, nmark, center_mark, nmark);
		
		LIBMVL_OFFSET64 j=mvl_randomize_bits64(nmark) & (N2-1);
		long long k;

		
		//fprintf(stderr, "%lld %lld\n", query_mark[i], j);
		
		k=first_mark[j];
		while(k>=0) {
			if(mark[k-1]==nmark)break;
			k=prev_mark[k-1];
			}
		if(k<0)continue;
		
	//	fprintf(stderr, "%lld %lld %lld\n", query_mark[i], j, k);
		
		for(long long m=first[k-1]; m>=0 ; m=prev[m-1]) {
			pd[indices_free]=m;
			indices_free++;
			}
			
		if(indices_free>indices_size) {
			Rprintf("*** INTERNAL ERROR: array overflow");
			}
		}

	//SETLENGTH(sa, indices_free);
	
	sa2=PROTECT(allocVector(REALSXP, indices_free));
	// ENABLE_REFCNT(sa);
	// INCREMENT_REFCNT(sa);
	// INCREMENT_NAMED(sa);
	pd2=REAL(sa2);
	
	for(LIBMVL_OFFSET64 m=0;m<indices_free;m++) pd2[m]=pd[m];
	
	SETCADR(R_fcall, ScalarReal(i+1));
	SETCADDR(R_fcall, sa2);
//	fprintf(stderr, "%lld %d %d %d (a)\n", i, MAYBE_REFERENCED(sa), -1, -1);
	tmp=PROTECT(eval(R_fcall, env));
//	if(MAYBE_REFERENCED(tmp))tmp=duplicate(tmp);
//	fprintf(stderr, "%lld %d %d %d (b)\n", i, MAYBE_REFERENCED(sa), -1, MAYBE_REFERENCED(tmp));

	SET_VECTOR_ELT(ans, i, tmp);
	UNPROTECT(2);
	}
	
free(values);
free(query_mark);
free(indices);
free(ball);
mvl_free_named_list(L);

UNPROTECT(3);
return(ans);
}


SEXP get_groups(SEXP prev, SEXP indices) 
{
LIBMVL_OFFSET64 count, data_offset, k, N, *v_idx, Nv;
int data_idx;
LIBMVL_VECTOR *vec;
SEXP ans;
double *apd;

decode_mvl_object(prev, &data_idx, &data_offset);
vec=get_mvl_vector(data_idx, data_offset);

if(vec==NULL) {
	error("Not an MVL object");
	return(R_NilValue);
	}
	
Nv=mvl_vector_length(vec);
	
if(get_indices(indices, vec, &N, &v_idx)) {
	return(R_NilValue);		
	}
	
count=0;
switch(mvl_vector_type(vec)) {
	case LIBMVL_VECTOR_INT64: {
		long long *pi=mvl_vector_data_int64(vec);
		for(LIBMVL_OFFSET64 i=0;i<N;i++) {
			k=v_idx[i];
			while(k>=0 && k<Nv) {
				k=pi[k]-1;
				count++;
				}
			}
		break;
		}
	case LIBMVL_VECTOR_INT32: {
		int *pi=mvl_vector_data_int32(vec);
		for(LIBMVL_OFFSET64 i=0;i<N;i++) {
			k=v_idx[i];
			while(k>=0 && k<Nv) {
				k=pi[k]-1;
				count++;
				}
			}
		break;
		}
	case LIBMVL_VECTOR_DOUBLE: {
		double *pi=mvl_vector_data_double(vec);
		for(LIBMVL_OFFSET64 i=0;i<N;i++) {
			k=v_idx[i];
			while(k>=0 && k<Nv) {
				k=pi[k]-1;
				count++;
				}
			}
		break;
		}
	default:
		error("Cannot process MVL vector of type %d\n", mvl_vector_type(vec));
		free(v_idx);
		return(R_NilValue);
	}
	
ans=PROTECT(allocVector(REALSXP, count));
apd=REAL(ans);

count=0;
switch(mvl_vector_type(vec)) {
	case LIBMVL_VECTOR_INT64: {
		long long *pi=mvl_vector_data_int64(vec);
		for(LIBMVL_OFFSET64 i=0;i<N;i++) {
			k=v_idx[i];
			while(k>=0 && k<Nv) {
				apd[count]=k+1;
				k=pi[k]-1;
				count++;
				}
			}
		break;
		}
	case LIBMVL_VECTOR_INT32: {
		int *pi=mvl_vector_data_int32(vec);
		for(LIBMVL_OFFSET64 i=0;i<N;i++) {
			k=v_idx[i];
			while(k>=0 && k<Nv) {
				apd[count]=k+1;
				k=pi[k]-1;
				count++;
				}
			}
		break;
		}
	case LIBMVL_VECTOR_DOUBLE: {
		double *pi=mvl_vector_data_double(vec);
		for(LIBMVL_OFFSET64 i=0;i<N;i++) {
			k=v_idx[i];
			while(k>=0 && k<Nv) {
				apd[count]=k+1;
				k=pi[k]-1;
				count++;
				}
			}
		break;
		}
	default:
		error("Cannot process MVL vector of type %d\n", mvl_vector_type(vec));
		free(v_idx);
		return(R_NilValue);
	}

free(v_idx);
UNPROTECT(1);
return(ans);
}

SEXP find_matches(SEXP data_list0, SEXP indices0, SEXP data_list1, SEXP indices1)
{
int data_idx;

LIBMVL_OFFSET64 data_offset, i, pairs_size;

double *pd;
int err;

void **vec_data0, **vec_data1;
LIBMVL_VECTOR **vectors0, **vectors1;
LIBMVL_OFFSET64 *v_idx0, *v_idx1, *key_hash, *key_last, *key_match_indices, *match_indices, *vec_data0_length, *vec_data1_length;
LIBMVL_OFFSET64 N0, N1;

HASH_MAP *hm;

SEXP ans, obj;
	
if(TYPEOF(data_list0)!=VECSXP) {
	error("order_vectors first argument must be a list of data to merge");
	return(R_NilValue);
	}

if(TYPEOF(data_list1)!=VECSXP) {
	error("order_vectors third argument must be a list of data to merge");
	return(R_NilValue);
	}

if(xlength(data_list0)<1 || xlength(data_list1)<1) {
	error("Vector lists should not be empty");
	return(R_NilValue);
	}

if(xlength(data_list0)!=xlength(data_list1)) {
	error("Vector lists should have the same number of vectors");
	return(R_NilValue);
	}
	
if(TYPEOF(indices0)!=NILSXP && xlength(indices0)<1) {
	error("Nothing to merge");
	return(R_NilValue);
	}

if(TYPEOF(indices1)!=NILSXP && xlength(indices1)<1) {
	error("Nothing to merge");
	return(R_NilValue);
	}
	
//Rprintf("Allocating vectors\n");
	
vec_data0=calloc(xlength(data_list0), sizeof(*vec_data0));
vec_data0_length=calloc(xlength(data_list0), sizeof(*vec_data0_length));
vectors0=calloc(xlength(data_list0), sizeof(*vectors0));
vec_data1=calloc(xlength(data_list1), sizeof(*vec_data1));
vec_data1_length=calloc(xlength(data_list1), sizeof(*vec_data1_length));
vectors1=calloc(xlength(data_list1), sizeof(*vectors1));
if(vec_data0==NULL || vec_data0_length==NULL || vectors0==NULL || vec_data1==NULL || vec_data1_length==NULL || vectors1==NULL) {
	error("Not enough memory");
	free(vec_data0);
	free(vectors0);
	free(vec_data1);
	free(vectors1);
	free(vec_data0_length);
	free(vec_data1_length);
	return(R_NilValue);
	}
	
//Rprintf("Computing data lists\n");
for(LIBMVL_OFFSET64 k=0;k<xlength(data_list0);k++) {
	decode_mvl_object(PROTECT(VECTOR_ELT(data_list0, k)), &data_idx, &data_offset);
	UNPROTECT(1);
	vectors0[k]=get_mvl_vector(data_idx, data_offset);
	
	if(vectors0[k]==NULL) {
		error("Invalid MVL object in first data list");
		free(vec_data0);
		free(vectors0);
		free(vec_data1);
		free(vectors1);
		free(vec_data0_length);
		free(vec_data1_length);
		return(R_NilValue);
		}
		
	vec_data0[k]=libraries[data_idx].data;
	vec_data0_length[k]=libraries[data_idx].length;
	
	decode_mvl_object(PROTECT(VECTOR_ELT(data_list1, k)), &data_idx, &data_offset);
	UNPROTECT(1);
	vectors1[k]=get_mvl_vector(data_idx, data_offset);
	
	if(vectors1[k]==NULL) {
		error("Invalid MVL object in second data list");
		free(vec_data0);
		free(vectors0);
		free(vec_data1);
		free(vectors1);
		free(vec_data0_length);
		free(vec_data1_length);
		return(R_NilValue);
		}
		
	vec_data1[k]=libraries[data_idx].data;
	vec_data1_length[k]=libraries[data_idx].length;
	
	if(mvl_vector_type(vectors0[k])!=mvl_vector_type(vectors1[k])) {
		error("Vector types do not match");
		free(vec_data0);
		free(vectors0);
		free(vec_data1);
		free(vectors1);
		free(vec_data0_length);
		free(vec_data1_length);
		return(R_NilValue);
		}
	}
	
//Rprintf("Extracting index0\n");
if(get_indices(indices0, vectors0[0], &N0, &v_idx0)) {
	free(vec_data0);
	free(vectors0);
	free(vec_data1);
	free(vectors1);
	free(vec_data0_length);
	free(vec_data1_length);
	return(R_NilValue);		
	}

//Rprintf("Extracting index1\n");
if(get_indices(indices1, vectors1[0], &N1, &v_idx1)) {
	free(vec_data0);
	free(vectors0);
	free(vec_data1);
	free(vectors1);
	free(v_idx0);
	free(vec_data0_length);
	free(vec_data1_length);
	return(R_NilValue);		
	}
	
//Rprintf("Computing key hash\n");

key_hash=calloc(N0, sizeof(*key_hash));
if(key_hash==NULL) {
	error("Not enough memory");
	free(vec_data0);
	free(vectors0);
	free(vec_data1);
	free(vectors1);
	free(v_idx0);
	free(v_idx1);
	free(vec_data0_length);
	free(vec_data1_length);
	return(R_NilValue);
	}

if((err=mvl_hash_indices(N0, v_idx0, key_hash, xlength(data_list0), vectors0, vec_data0, vec_data0_length, LIBMVL_COMPLETE_HASH))!=0) {
	free(vec_data0);
	free(vectors0);
	free(vec_data1);
	free(vectors1);
	free(v_idx0);
	free(v_idx1);
	free(key_hash);
	free(vec_data0_length);
	free(vec_data1_length);
	error("Error hashing key indices %d\n", err);
	return(R_NilValue);
	}
	
// Rprintf("Allocating hash map\n");
hm=mvl_allocate_hash_map(N1);
hm->hash_count=N1;

//Rprintf("Computing data hash\n");
if((err=mvl_hash_indices(N1, v_idx1, hm->hash, xlength(data_list1), vectors1, vec_data1, vec_data1_length, LIBMVL_COMPLETE_HASH))!=0) {
	free(vec_data0);
	free(vectors0);
	free(vec_data1);
	free(vectors1);
	free(v_idx0);
	free(v_idx1);
	free(key_hash);
	free(vec_data0_length);
	free(vec_data1_length);
	mvl_free_hash_map(hm);
	error("Error hashing indices %d\n", err);
	return(R_NilValue);
	}
	
//Rprintf("Computing hash map\n");
mvl_compute_hash_map(hm);

//Rprintf("Estimating match count\n");
pairs_size=mvl_hash_match_count(N0, key_hash, hm);

if(pairs_size>1e9) {
	Rprintf("Expecting %lld matches\n", pairs_size);
	}
	
key_last=calloc(N0, sizeof(*key_last));
key_match_indices=calloc(pairs_size, sizeof(*key_match_indices));
match_indices=calloc(pairs_size, sizeof(*match_indices));

if(key_last==NULL || key_match_indices==NULL || match_indices==NULL) {
	free(vec_data0);
	free(vectors0);
	free(vec_data1);
	free(vectors1);
	free(v_idx0);
	free(v_idx1);
	free(key_hash);
	free(key_last);
	free(key_match_indices);
	free(match_indices);
	free(vec_data0_length);
	free(vec_data1_length);
	mvl_free_hash_map(hm);
	error("Not enough memory");
	return(R_NilValue);
	}

//Rprintf("Finding matches\n");
if((err=mvl_find_matches(N0, v_idx0, xlength(data_list0), vectors0, vec_data0, vec_data0_length, key_hash,
	N1, v_idx1, xlength(data_list1), vectors1, vec_data1, vec_data1_length, hm,
	key_last, pairs_size, key_match_indices, match_indices))) {
	error("Error computing merge plan %d\n", err);
	free(vec_data0);
	free(vectors0);
	free(vec_data1);
	free(vectors1);
	free(v_idx0);
	free(v_idx1);
	free(key_hash);
	free(key_last);
	free(key_match_indices);
	free(match_indices);
	free(vec_data0_length);
	free(vec_data1_length);
	mvl_free_hash_map(hm);
	return(R_NilValue);
	}

// Rprintf("Formating results\n");
	
mvl_free_hash_map(hm);
free(key_hash);
	
ans=PROTECT(allocVector(VECSXP, 3));

obj=PROTECT(allocVector(REALSXP, N0+1));
pd=REAL(obj);
pd[0]=1;

for(i=0;i<N0;i++) pd[i+1]=key_last[i]+1;

SET_VECTOR_ELT(ans, 0, obj);
UNPROTECT(1);

obj=PROTECT(allocVector(REALSXP, key_last[N0-1]));
pd=REAL(obj);

for(i=0;i<key_last[N0-1];i++) pd[i]=key_match_indices[i]+1;

SET_VECTOR_ELT(ans, 1, obj);
UNPROTECT(1);

obj=PROTECT(allocVector(REALSXP, key_last[N0-1]));
pd=REAL(obj);

for(i=0;i<key_last[N0-1];i++) pd[i]=match_indices[i]+1;

SET_VECTOR_ELT(ans, 2, obj);
UNPROTECT(1);
UNPROTECT(1);

free(vec_data0);
free(vectors0);
free(vec_data1);
free(vectors1);
free(vec_data0_length);
free(vec_data1_length);
free(v_idx0);
free(v_idx1);
free(key_last);
free(key_match_indices);
free(match_indices);
return(ans);
}

SEXP mvl_xlength_int(SEXP obj)
{
SEXP ans;

ans=PROTECT(allocVector(REALSXP, 1));
REAL(ans)[0]=xlength(obj);
UNPROTECT(1);
return(ans);
}

SEXP group_vectors(SEXP data_list, SEXP indices)
{
SEXP ans, obj, data;
int err;
double *pd, *fd;

void **vec_data;
LIBMVL_VECTOR **vectors;
LIBMVL_OFFSET64 *v_idx, *vec_data_length;
LIBMVL_OFFSET64 N;
int data_idx;

LIBMVL_OFFSET64 data_offset, i, j, k;

HASH_MAP *hm;

	
if(TYPEOF(data_list)!=VECSXP) {
	error("group_vectors first argument must be a list of data to group");
	return(R_NilValue);
	}


if(xlength(data_list)<1) {
	error("Vector lists should not be empty");
	return(R_NilValue);
	}

if(TYPEOF(indices)!=NILSXP && xlength(indices)<1) {
	error("Nothing to group");
	return(R_NilValue);
	}

//Rprintf("Allocating vectors\n");
	
vec_data=calloc(xlength(data_list), sizeof(*vec_data));
vec_data_length=calloc(xlength(data_list), sizeof(*vec_data_length));
vectors=calloc(xlength(data_list), sizeof(*vectors));
if(vec_data==NULL || vec_data_length==NULL || vectors==NULL) {
	error("Not enough memory");
	free(vec_data);
	free(vec_data_length);
	free(vectors);
	return(R_NilValue);
	}
	
//Rprintf("Computing data lists\n");
for(LIBMVL_OFFSET64 k=0;k<xlength(data_list);k++) {
	decode_mvl_object(PROTECT(VECTOR_ELT(data_list, k)), &data_idx, &data_offset);
	UNPROTECT(1);
	vectors[k]=get_mvl_vector(data_idx, data_offset);
	
	if(vectors[k]==NULL) {
		error("Invalid MVL object in first data list");
		free(vec_data);
		free(vec_data_length);
		free(vectors);
		return(R_NilValue);
		}
		
	vec_data[k]=libraries[data_idx].data;
	vec_data_length[k]=libraries[data_idx].length;
	}
	
//Rprintf("Extracting index\n");
if(get_indices(indices, vectors[0], &N, &v_idx)) {
	free(vec_data);
	free(vectors);
	return(R_NilValue);		
	}
	
hm=mvl_allocate_hash_map(N);
hm->hash_count=N;

//Rprintf("Computing data hash\n");
if((err=mvl_hash_indices(N, v_idx, hm->hash, xlength(data_list), vectors, vec_data, vec_data_length, LIBMVL_COMPLETE_HASH))!=0) {
	free(vec_data);
	free(vec_data_length);
	free(vectors);
	free(v_idx);
	mvl_free_hash_map(hm);
	error("Error hashing indices %d\n", err);
	return(R_NilValue);
	}
	
//Rprintf("Computing hash map hash_size=%lld hash_map_size=%lld\n", hm->hash_size, hm->hash_map_size);
mvl_compute_hash_map(hm);

//Rprintf("Finding groups\n");
mvl_find_groups(N, v_idx, xlength(data_list), vectors, vec_data, vec_data_length, hm);

ans=PROTECT(allocVector(VECSXP, 2));

obj=PROTECT(allocVector(REALSXP, N));
data=PROTECT(allocVector(REALSXP, hm->first_count+1));

pd=REAL(obj);
fd=REAL(data);
fd[0]=1;

j=0;
for(i=0;i<hm->first_count;i++) {
	k=hm->first[i];
	while(k!=~0LLU) {
		pd[j]=v_idx[k]+1;
		j++;
		k=hm->next[k];
		}
	fd[i+1]=j+1;
	}
	
SET_VECTOR_ELT(ans, 0, data);
SET_VECTOR_ELT(ans, 1, obj);

free(vec_data);
free(vec_data_length);
free(vectors);
free(v_idx);
mvl_free_hash_map(hm);

UNPROTECT(3);
return(ans);
}

SEXP group_lapply(SEXP si, SEXP index, SEXP fn, SEXP env)
{
LIBMVL_OFFSET64 N, Nidx, stretch_max;
SEXP ans, vidx, R_fcall;

double *psi, *pidx, *pidx2;

if(xlength(si)<2) {
	error("stretch index should have length of at least 2");
	return(R_NilValue);
	}
if(!isFunction(fn)) {
	error("third argument must be a function");
	return(R_NilValue);
	}
if(!isEnvironment(env)) {
	error("fourth argument should be an environment");
	return(R_NilValue);
	}
	
N=xlength(si)-1;
psi=REAL(si);

Nidx=xlength(index);
pidx=REAL(index);

ans=PROTECT(allocVector(VECSXP, N));

R_fcall = PROTECT(lang2(fn, R_NilValue));

#if 0
stretch_max=1;
for(LIBMVL_OFFSET64 i=0;i<N;i++) {
	LIBMVL_OFFSET64 k=psi[i+1]-psi[i];
	if(k>stretch_max)stretch_max=k;
	}
	
vidx=PROTECT(allocVector(REALSXP, stretch_max));
pidx2=REAL(vidx);
#endif

for(LIBMVL_OFFSET64 i=0;i<N;i++) {
	LIBMVL_OFFSET64 k0, k1;
	k0=psi[i]-1;
	k1=psi[i+1]-1;
	if(k1<=k0)continue;
	if(k0 >= Nidx || k1>Nidx)continue;
	//SETLENGTH(vidx, k1-k0);

	vidx=PROTECT(allocVector(REALSXP, k1-k0));
	pidx2=REAL(vidx);
	
	for(LIBMVL_OFFSET64 j=k0;j<k1;j++) {
		pidx2[j-k0]=pidx[j];
		}
	SETCADR(R_fcall, vidx);
	SET_VECTOR_ELT(ans, i, PROTECT(eval(R_fcall, env)));
	UNPROTECT(2);
	}

UNPROTECT(2);
return(ans);
}


SEXP compute_repeats(SEXP data_list)
{
int data_idx;
LIBMVL_OFFSET64 data_offset;

void **vec_data;
LIBMVL_OFFSET64 *vec_data_length;
LIBMVL_VECTOR **vectors;
double *dans;

LIBMVL_PARTITION el;

SEXP ans;
	

if(TYPEOF(data_list)!=VECSXP) {
	error("compute_repeats first argument must be a list of data to sort");
	return(R_NilValue);
	}

if(xlength(data_list)<1) {
	error("No hashes to compute");
	return(R_NilValue);
	}
	
vec_data=calloc(xlength(data_list), sizeof(*vec_data));
vec_data_length=calloc(xlength(data_list), sizeof(*vec_data_length));
vectors=calloc(xlength(data_list), sizeof(*vectors));
if(vec_data==NULL || vec_data_length==NULL || vectors==NULL) {
	error("Not enough memory");
	free(vec_data);
	free(vec_data_length);
	free(vectors);
	return(R_NilValue);
	}

for(LIBMVL_OFFSET64 k=0;k<xlength(data_list);k++) {
	decode_mvl_object(PROTECT(VECTOR_ELT(data_list, k)), &data_idx, &data_offset);
	UNPROTECT(1);
	vectors[k]=get_mvl_vector(data_idx, data_offset);
	
	if(vectors[k]==NULL) {
		error("Invalid MVL object in data list");
		free(vec_data);
		free(vec_data_length);
		free(vectors);
		return(R_NilValue);
		}
		
	vec_data[k]=libraries[data_idx].data;
	vec_data_length[k]=libraries[data_idx].length;
	}

memset(&el, 0, sizeof(el));
mvl_find_repeats(&el, xlength(data_list), vectors, vec_data, vec_data_length);

ans=PROTECT(allocVector(REALSXP, el.count+1));
dans=REAL(ans);

for(LIBMVL_OFFSET64 i=0;i<el.count;i++) dans[i]=el.offset[i]+1;

mvl_free_partition_arrays(&el);
free(vec_data);
free(vec_data_length);
free(vectors);
	
UNPROTECT(1);

return(ans);	
}

SEXP write_extent_index(SEXP idx0, SEXP data_list)
{
int idx, data_idx;
LIBMVL_OFFSET64 data_offset;

void **vec_data;
LIBMVL_VECTOR **vectors;
LIBMVL_OFFSET64 offset, *vec_data_length;
double *doffset=(double *)&offset;

LIBMVL_EXTENT_INDEX ei;

SEXP ans, class;
	
if(length(idx0)!=1) {
	error("write_extent_index first argument must be a single integer");
	return(R_NilValue);
	}
idx=INTEGER(idx0)[0];
if(idx<0 || idx>=libraries_free) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].ctx==NULL) {
	error("invalid MVL handle");
	return(R_NilValue);
	}
if(libraries[idx].f==NULL) {
	error("library not open for writing");
	return(R_NilValue);
	}
if(mvl_get_error(libraries[idx].ctx)!=0) {
	error("library has error status %d: %s", mvl_get_error(libraries[idx].ctx), mvl_strerror(libraries[idx].ctx));
	return(R_NilValue);
	}


if(TYPEOF(data_list)!=VECSXP) {
	error("compute_extent_index second argument must be a list of data to index");
	return(R_NilValue);
	}

if(xlength(data_list)<1) {
	error("No vectors to index");
	return(R_NilValue);
	}
	
vec_data=calloc(xlength(data_list), sizeof(*vec_data));
vec_data_length=calloc(xlength(data_list), sizeof(*vec_data_length));
vectors=calloc(xlength(data_list), sizeof(*vectors));
if(vec_data==NULL || vectors==NULL) {
	error("Not enough memory");
	free(vec_data);
	free(vec_data_length);
	free(vectors);
	return(R_NilValue);
	}

for(LIBMVL_OFFSET64 k=0;k<xlength(data_list);k++) {
	decode_mvl_object(PROTECT(VECTOR_ELT(data_list, k)), &data_idx, &data_offset);
	UNPROTECT(1);
	vectors[k]=get_mvl_vector(data_idx, data_offset);
	
	if(vectors[k]==NULL) {
		error("Invalid MVL object in data list");
		free(vec_data);
		free(vec_data_length);
		free(vectors);
		return(R_NilValue);
		}
		
	vec_data[k]=libraries[data_idx].data;
	vec_data_length[k]=libraries[data_idx].length;
	}

mvl_init_extent_index(&ei);
mvl_compute_extent_index(&ei, xlength(data_list), vectors, vec_data, vec_data_length);
offset=mvl_write_extent_index(libraries[idx].ctx, &ei);

mvl_free_extent_index_arrays(&ei);

free(vec_data);
free(vec_data_length);
free(vectors);
	
ans=PROTECT(allocVector(REALSXP, 1));
REAL(ans)[0]=*doffset;

class=PROTECT(allocVector(STRSXP, 1));
SET_STRING_ELT(class, 0, mkChar("MVL_OFFSET"));
classgets(ans, class);
UNPROTECT(2);

return(ans);	
}

SEXP extent_index_lapply(SEXP extent_index0, SEXP data_list, SEXP fn, SEXP env)
{
LIBMVL_OFFSET64 data_offset, index_offset, Nv, indices_free, *hash, k;
int data_idx, index_idx, err;
LIBMVL_VECTOR *vec;
SEXP ans, sa, R_fcall, tmp;
double *pd;
LIBMVL_EXTENT_INDEX ei;
LIBMVL_EXTENT_LIST el;

if(TYPEOF(data_list)!=VECSXP) {
	error("Second argument should be a list (or data frame) of vectors to query");
	return(R_NilValue);
	}
	
decode_mvl_object(extent_index0, &index_idx, &index_offset);

if(index_idx<0) {
	error("extent index is not an MVL OBJECT");
	return(R_NilValue);
	}

mvl_init_extent_index(&ei);

if((err=mvl_load_extent_index(libraries[index_idx].ctx, libraries[index_idx].data, libraries[index_idx].length, index_offset, &ei))!=0) {
	error("Error accessing extent index (%d): %s", err, mvl_strerror(libraries[index_idx].ctx));
	return(R_NilValue);
	}

if(xlength(data_list)!=ei.hash_map.vec_count) {
	error("Number of vectors (columns) does not match - original index used %lld vectors", ei.hash_map.vec_count);
	mvl_free_extent_index_arrays(&ei);
	return(R_NilValue);
	}
	
switch(TYPEOF(VECTOR_ELT(data_list, 0))) {
	case REALSXP:
	case INTSXP:
	case STRSXP:
		Nv=xlength(VECTOR_ELT(data_list, 0));
		break;
	case VECSXP:
		decode_mvl_object(PROTECT(VECTOR_ELT(data_list, 0)), &data_idx, &data_offset);
		UNPROTECT(1);
		vec=get_mvl_vector(data_idx, data_offset);
		if(vec==NULL) {
			mvl_free_extent_list_arrays(&el);
			error("Not an MVL object");
			return(R_NilValue);
			}
		Nv=mvl_vector_length(vec);
		if(mvl_vector_type(vec)==LIBMVL_PACKED_LIST64)Nv--;
		break;
	default:
		mvl_free_extent_index_arrays(&ei);
		error("Cannot handle R vector of type %d", TYPEOF(VECTOR_ELT(data_list, 0)));
		return(R_NilValue);
		break;
	}
	

hash=calloc(Nv, sizeof(*hash));

for(LIBMVL_OFFSET64 i=0;i<Nv;i++)hash[i]=MVL_SEED_HASH_VALUE;

for(LIBMVL_OFFSET64 i=0;i<xlength(data_list);i++) {
	if(hash_vector_range(PROTECT(VECTOR_ELT(data_list, i)), 0, Nv, ei.hash_map.vec_types[i], hash)) {
		UNPROTECT(1);
		return(R_NilValue);
		}
	UNPROTECT(1);
	}
	
for(LIBMVL_OFFSET64 i=0;i<Nv;i++)hash[i]=mvl_randomize_bits64(hash[i]);
	
ans=PROTECT(allocVector(VECSXP, Nv));	

R_fcall = PROTECT(lang3(fn, R_NilValue, R_NilValue));
	
mvl_init_extent_list(&el);

for(LIBMVL_OFFSET64 i=0;i<Nv;i++) {
	mvl_empty_extent_list(&el);
	mvl_get_extents(&ei, hash[i], &el);
	
	indices_free=0;
	
	for(LIBMVL_OFFSET64 j=0;j<el.count;j++)indices_free+=el.stop[j]-el.start[j];
	
	if(indices_free<1)continue;
	
	sa=PROTECT(allocVector(REALSXP, indices_free));
	pd=REAL(sa);
	
	k=0;
	for(LIBMVL_OFFSET64 j=0;j<el.count;j++) {
		for(LIBMVL_OFFSET64 m=el.start[j]; m<el.stop[j];m++,k++)
			pd[k]=m+1;
		}
	

	SETCADR(R_fcall, ScalarReal(i+1));
	SETCADDR(R_fcall, sa);
//	fprintf(stderr, "%lld %d %d %d (a)\n", i, MAYBE_REFERENCED(sa), -1, -1);
	tmp=PROTECT(eval(R_fcall, env));
//	if(MAYBE_REFERENCED(tmp))tmp=duplicate(tmp);
//	fprintf(stderr, "%lld %d %d %d (b)\n", i, MAYBE_REFERENCED(sa), -1, MAYBE_REFERENCED(tmp));

	SET_VECTOR_ELT(ans, i, tmp);
	UNPROTECT(2);
	}
	
mvl_free_extent_list_arrays(&el);
	
free(hash);

UNPROTECT(2);
return(ans);
}

SEXP extent_index_scan(SEXP extent_index0, SEXP fn, SEXP env)
{
LIBMVL_OFFSET64 index_offset, Nv, indices_free, k;
int index_idx, err;
SEXP ans, sa, R_fcall, tmp;
double *pd;
LIBMVL_EXTENT_INDEX ei;
LIBMVL_EXTENT_LIST el;

decode_mvl_object(extent_index0, &index_idx, &index_offset);

mvl_init_extent_index(&ei);

if((err=mvl_load_extent_index(libraries[index_idx].ctx, libraries[index_idx].data, libraries[index_idx].length, index_offset, &ei))!=0) {
	error("Error accessing extent index (%d): %s", err, mvl_strerror(libraries[index_idx].ctx));
	return(R_NilValue);
	}


Nv=ei.hash_map.hash_count;

	
ans=PROTECT(allocVector(VECSXP, Nv));	

R_fcall = PROTECT(lang3(fn, R_NilValue, R_NilValue));
	
mvl_init_extent_list(&el);

for(LIBMVL_OFFSET64 i=0;i<Nv;i++) {
	mvl_empty_extent_list(&el);
	mvl_get_extents(&ei, ei.hash_map.hash[i], &el);
	
	indices_free=0;
	
	for(LIBMVL_OFFSET64 j=0;j<el.count;j++)indices_free+=el.stop[j]-el.start[j];
	
	if(indices_free<1)continue;
	
	sa=PROTECT(allocVector(REALSXP, indices_free));
	pd=REAL(sa);
	
	k=0;
	for(LIBMVL_OFFSET64 j=0;j<el.count;j++) {
		for(LIBMVL_OFFSET64 m=el.start[j]; m<el.stop[j];m++,k++)
			pd[k]=m+1;
		}
	

	SETCADR(R_fcall, ScalarReal(i+1));
	SETCADDR(R_fcall, sa);
//	fprintf(stderr, "%lld %d %d %d (a)\n", i, MAYBE_REFERENCED(sa), -1, -1);
	tmp=PROTECT(eval(R_fcall, env));
//	if(MAYBE_REFERENCED(tmp))tmp=duplicate(tmp);
//	fprintf(stderr, "%lld %d %d %d (b)\n", i, MAYBE_REFERENCED(sa), -1, MAYBE_REFERENCED(tmp));

	SET_VECTOR_ELT(ans, i, tmp);
	UNPROTECT(2);
	}
	
mvl_free_extent_list_arrays(&el);
mvl_free_extent_index_arrays(&ei);
	
UNPROTECT(2);
return(ans);
}


static const R_CMethodDef cMethods[] = {
//    {"foo", (DL_FUNC) &foo, 4, foo_t},
//    {"bar_sym", (DL_FUNC) &bar, 0},
   {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[] = {
   {"mmap_library", (DL_FUNC) &mmap_library, 2},
   {"remap_library", (DL_FUNC) &remap_library, 2},
  {"close_library",  (DL_FUNC) &close_library, 1},
  {"find_directory_entries",  (DL_FUNC) &find_directory_entries, 2},
  {"get_directory",  (DL_FUNC) &get_directory, 1},
  {"read_metadata",  (DL_FUNC) &read_metadata, 2},
  {"read_lengths",  (DL_FUNC) &read_lengths, 2},
  {"compute_vector_stats",  (DL_FUNC) &compute_vector_stats, 2},
  {"read_types",  (DL_FUNC) &read_types, 2},
  {"get_vector_data_ptr",  (DL_FUNC) &get_vector_data_ptr, 2},
  {"read_vectors_raw",  (DL_FUNC) &read_vectors_raw, 2},
  {"read_vectors_idx_raw2",  (DL_FUNC) &read_vectors_idx_raw2, 3},
  {"read_vectors",  (DL_FUNC) &read_vectors, 2},
  {"read_vectors_idx3",  (DL_FUNC) &read_vectors_idx3, 3},
  {"add_directory_entries",  (DL_FUNC) &add_directory_entries, 3},
  {"write_vector",  (DL_FUNC) &write_vector, 4},
  {"start_write_vector",  (DL_FUNC) &start_write_vector, 5},
  {"rewrite_vector",  (DL_FUNC) &rewrite_vector, 3},
  {"fused_write_vector",  (DL_FUNC) &fused_write_vector, 4},
  {"order_vectors",  (DL_FUNC) &order_vectors, 3},
  {"hash_vectors",  (DL_FUNC) &hash_vectors, 2},
  {"write_hash_vectors",  (DL_FUNC) &write_hash_vectors, 2},
  {"find_matches",  (DL_FUNC) &find_matches, 4},
  {"indexed_copy_vector",  (DL_FUNC) &indexed_copy_vector, 4},
  {"mvl_xlength_int",  (DL_FUNC) &mvl_xlength_int, 1},
  {"group_vectors",  (DL_FUNC) &group_vectors, 2},
  {"group_lapply",  (DL_FUNC) &group_lapply, 4},
  {"write_groups",  (DL_FUNC) &write_groups, 2},
  {"get_groups",  (DL_FUNC) &get_groups, 2},
  {"write_spatial_groups",  (DL_FUNC) &write_spatial_groups, 3},
  {"get_neighbors",  (DL_FUNC) &get_neighbors, 2},
  {"neighbors_lapply",  (DL_FUNC) &neighbors_lapply, 4},
  {"write_extent_index",  (DL_FUNC) &write_extent_index, 2},
  {"compute_repeats",  (DL_FUNC) &compute_repeats, 1},
  {"extent_index_lapply",  (DL_FUNC) &extent_index_lapply, 4},
  {"extent_index_scan",  (DL_FUNC) &extent_index_scan, 3},
  {"get_status",  (DL_FUNC) &get_status, 0},
   {NULL, NULL, 0}
};



void R_init_RMVL(DllInfo *info) {

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);

}
