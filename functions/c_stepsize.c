#include "mex.h"
#include "math.h"
#include "time.h"
#include "stdlib.h"

/* get alpha by sweeping,  minimise median residual of cameras */

double eps_feas = 0.000000000000001; /* 1e-15 */
double eps_root = 0.0000000001; /* 1e-10*/
double eps_diff = 0.00000001; /* 1e-8 */
double Inf = 100000000;
int compare(const void *a,const void *b);
int compare_descent(const void *a,const void *b);
double* getbound(double* a, double* b, double* c, double* d, int M, double* X, double* dir);
void getorder(double* a, double* b, double* c, double* d, int M, double* X, int K, int* List, int* aid, int* anchor);
double alpha_bf(double Ps[][3], int i);
double alpha_af(double Ps[][3], int i, int M);

int tangent(double* a, double* b, double* c, double* d, \
            int M, double* X, double* dir, \
            double a_bf, double a_af, int cam1, int cam2);

double err_of_root(double* a1, double b1, double* c1, double d1, double* X, double* dir, double rt);

/* get error of camera cam for given list of alphas */
void get_error(double* a, double* b, double* c, double* d, \
               double* X, double* dir, \
               int M, double* alphas, double* errors, int n_alf, int cam);

/* get extremity */
void get_extreme(double* a, double* b, double* c, double* d, \
                   double lb, double ub, double* X, double* dir, \
                   int M, double Ps[][3], int* Ps_i);
   
/* get intersections */
void get_intersection(double* a, double* b, double* c, double* d, \
                   double lb, double ub, double* X, double* dir, \
                   int M, double Ps[][3], int* Ps_i);

/* get positive roots of constraints */
void f_root_positive( double* a, double* b, double* c, double* d, \
                        int cid_i, int* cid_j, int id_length, \
                        double* X, double* dir, double lb, double ub, \
                        int M, double *roots, double* vs, int* n_rt );



/* double* roots getroots(double** aa, double* bb, double** cc, double* dd, int i, int j, double *X, double* dir, double K);*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int debug = 0;
    double *a, *b, *c, *d, *e, *f, *g;
    double *alpha;
    int m1,n1,m2,n2,m3,n3,m4,n4,m5,n5,m6,n6;
    int M;
    
    /* Check for proper number of arguments. */    
    if (nrhs != 7)
    {
        mexErrMsgTxt("7 inputs required.");
    }
    else if (nlhs > 3)
    {
        mexErrMsgTxt("Too many output arguments.");
    }
  
    /* Assign pointers to inputs.*/
    a = mxGetPr(prhs[0]);
    b = mxGetPr(prhs[1]);
    c = mxGetPr(prhs[2]);
    d = mxGetPr(prhs[3]);
    e = mxGetPr(prhs[4]);
    f = mxGetPr(prhs[5]);
    g = mxGetPr(prhs[6]);
  
    /* Get sizes of input matrices.*/
    m1 = mxGetM(prhs[0]);
    n1 = mxGetN(prhs[0]);
    m2 = mxGetM(prhs[1]);
    n2 = mxGetN(prhs[1]);
    m3 = mxGetM(prhs[2]);
    n3 = mxGetN(prhs[2]);
    m4 = mxGetM(prhs[3]);
    n4 = mxGetN(prhs[3]);
    m5 = mxGetM(prhs[4]);
    n5 = mxGetN(prhs[4]);
    m6 = mxGetM(prhs[5]);
    n6 = mxGetN(prhs[5]);

 
    
    M = m1; /* number of constraints */
    int nCam = ceil(M/4);
    int K = (int)(*g); /* K normally is the index of the median */
    double X[3] = {e[0], e[1], e[2]};
    

    double dir[3] = {f[0], f[1], f[2]};
    
    int i,j,t;
    if(debug)
    {
        printf("X\t\tdir\t\t\n");
        for(t=0;t<3;t++)
        {
            printf("%f\t\t%f\t\t\n",X[t],dir[t]);
        }
    }    

    /* Create matrices for the return argument. */    
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(3*M*(M-1),1,mxREAL); 

    /* Assign pointers to output.*/
    /*Ps = mxGetPr(plhs[0]);*/
    alpha = mxGetPr(plhs[0]);
    double *Ps2 = mxGetPr(plhs[1]);
    
    
    /*-----------       get roots    ---------------*/
    double a1[3]; 
    double b1;    
    double c1[3]; 
    double d1;    
    double a2[3]; 
    double b2;    
    double c2[3]; 
    double d2;    

    double A,B,C;
    int Ps_i = 0;
    double temp;
    double rt;
    double Ps[M*(M-1)][3]; /* each row: [alpha, line_1_aid, line_2_aid]*/
    double *bound;
    double lb, ub;
    
    
    /* get bound;*/
    bound = getbound(a, b, c, d, M, X, dir);
    lb = bound[0];
    ub = bound[1];
    free(bound);
    
   
    /* get extremity */ 
    get_extreme(a, b, c, d, lb, ub, X, dir, M, Ps, &Ps_i);

    /* get intersections  */
    get_intersection(a, b, c, d, lb, ub, X, dir, M, Ps, &Ps_i);
    
    for(i=0;i<Ps_i;i++)
    {
        Ps2[3*i+0] = Ps[i][0];
        Ps2[3*i+1] = Ps[i][1];
        Ps2[3*i+2] = Ps[i][2];
    }
    /* end of get roots */


    
    qsort(Ps,Ps_i,3*sizeof(double),compare);

    for(i=0;i<Ps_i;i++)
    {
        Ps2[3*i+0] = Ps[i][0];
        Ps2[3*i+1] = Ps[i][1];
        Ps2[3*i+2] = Ps[i][2];
    }

    
    /*------------       Initialise List     ----------*/
    
    int List[nCam];  /*/ [order]*/
    int aid[nCam];
    int anchor; /*/ anchor == the K_th line_id*/
    double alf_s = Ps[0][0]/2;
    double Xs[3] = {X[0]+alf_s*dir[0], X[1]+alf_s*dir[1], X[2]+alf_s*dir[2]};
    
    getorder(a, b, c, d, M, Xs, K, List, aid, &anchor);

 
    /******************************  sweep  ******************************/
    double bestalpha = 0;
    double bestres = Inf;

    int id1, id2; /*/ list id (not the line id)*/
    double a_bf, a_af;
    double X1[3];
    double X_bf[3];
    double X_af[3];
    double res_bf[2];
    double res_af[2];
    int other;
    double inter; /*/ intermedia residual*/
    
    for(i=0; i<Ps_i; i++)
    {
        if(Ps[i][2]!=0) /* an intersection */
        {
            id1 = aid[(int)Ps[i][1]-1] - 1;
            id2 = aid[(int)Ps[i][2]-1] - 1;
            a_bf = alpha_bf(Ps, i);
            a_af = alpha_af(Ps, i, Ps_i);

            if(debug)
            {
                for(j=0;j<M;j++)
                {
                    printf("%d\t%d\t\n",aid[j],List[j]);
                }
                printf("anchor=%d\n\n",anchor);
                printf("%d\t%d\t%d\t%d\t\n",id1+1,id2+1,List[id1],List[id2]);
            }


            /*********************    check tangency   ********************/
            /*/ check if the two curves indeed cross*/
            if ( tangent(a, b, c, d, M, X, dir, a_bf, a_af, (int)Ps[i][1], (int)Ps[i][2]) )
            {
                /* printf("Tangent!!!  [%d, %d] at %f\n\n\n", (int)Ps[i][1], (int)Ps[i][2], Ps[i][0]); */
                continue;
            }

            /******************************  swap  ******************************/
            int temp = List[id1];
            List[id1] = List[id2];
            List[id2] = temp;
            temp = aid[(int)Ps[i][1]-1];
            aid[(int)Ps[i][1]-1] = aid[(int)Ps[i][2]-1];
            aid[(int)Ps[i][2]-1] = temp;
            if(debug)
            {
                printf("%d\t%d\t%d\t%d\t\n",id1+1,id2+1,List[id1],List[id2]);
                printf("\n");
            }
        }
        /******************************  check optimal  ******************************/
        
        if( anchor==(int)Ps[i][1] || anchor==(int)Ps[i][2] )
        {
            
            /**** firstly, update anchor ****/
            anchor = List[K-1];
        
            /**** secondly, check optimal ****/
            double alf_t = Ps[i][0];
            get_error(a, b, c, d, X, dir, M, &alf_t, &inter, 1, (int)Ps[i][1]);
            if(inter<=0)
                continue;
            if(inter < bestres)
            {
                bestalpha = Ps[i][0];
                bestres = inter;
            }
        }
        /****************************** END check optimal  ******************************/
    } /* end of for (Ps) */
    
    /****  need to check optimality for the last one intersection ****/
    X1[0] = X[0] + Ps[i-1][0]*dir[0]; 
    X1[1] = X[1] + Ps[i-1][0]*dir[1]; 
    X1[2] = X[2] + Ps[i-1][0]*dir[2]; 
    
    getorder(a, b, c, d, M, X1, K, List, aid, &anchor);  
    double alf_t = Ps[i-1][0];
    get_error(a, b, c, d, X, dir, M, &alf_t, &inter, 1, anchor);
    /*printf("inter = %.16f, bestres = %.16f\n", inter, bestres);*/
    if(inter < bestres)
    {
        /*printf("i = %d, updated with the last root %.16f\n",i, Ps[i][0]);*/
        bestalpha = Ps[i-1][0];
        bestres = inter;
    }  
    *alpha = bestalpha;
}


/* compare */
int compare(const void *a,const void *b)
{
    double *da = (double*) a;
    double *db = (double*) b;
    if (da[0] > db[0])
        return 1;
    else if(da[0] < db[0])
        return -1;
    else
    {        
        return 0;
          
        if (da[1] > db[1])
            return 1;
        else if(da[1] < db[1])
            return -1;
        else
        {
            if (da[2] > db[2])
                return 1;
            else if(da[2] < db[2])
                return -1;
            else
                return 0;
        }
    }
}

int compare_descent(const void *a,const void *b)
{
    double *da = (double*) a;
    double *db = (double*) b;
    if (db[0] > da[0])
        return 1;
    else if(db[0] < da[0])
        return -1;
    else
    {
        return 0;
                
        if (db[1] > da[1])
            return 1;
        else if(db[1] < da[1])
            return -1;
        else
        {
            if (db[2] > da[2])
                return 1;
            else if(db[2] < da[2])
                return -1;
            else
                return 0;
        }
    }
}

/* get the lower bound and the upper bound of feasible alpha*/
double* getbound(double* a, double* b, double* c, double* d, int M, double* X, double* dir)
{
    double lb = 0;
    double ub = Inf; /*/ just a large number;*/
    int i;
    double temp;
    double *result;
    double a1[3]; 
    double b1;    
    double c1[3]; 
    double d1;   
    
    for(i=0;i<M;i++)
    {
        a1[0] = a[i+0*M],  a1[1] = a[i+1*M],  a1[2] = a[i+2*M];
        c1[0] = c[i+0*M],  c1[1] = c[i+1*M],  c1[2] = c[i+2*M];
        b1 = b[i];
        d1 = d[i];
        temp = -(c1[0]*X[0] + c1[1]*X[1] + c1[2]*X[2] + d1) / (c1[0]*dir[0] + c1[1]*dir[1] + c1[2]*dir[2]);
        
        if( (c1[0]*dir[0] + c1[1]*dir[1] + c1[2]*dir[2]) > 0 ) /*/ for lower bound*/
        {
            if(temp > lb)
                lb = temp;
        }
        else /* for upper bound*/
        {
            if(temp < ub)
                ub = temp;
        }
    }
    result = (double*)malloc(2*sizeof(double));
    result[0] = lb;
    result[1] = ub;
    return result;
}


/* order curves into List and aid */
void getorder(double* a, double* b, double* c, double* d, int M, double* X, int K, int* List, int* aid, int* anchor)
{
    double a1[3]; 
    double b1;    
    double c1[3]; 
    double d1;    
    double a2[3]; 
    double b2;    
    double c2[3]; 
    double d2;
    int nCam = M/4;
    double reslist[M][2]; /* [error, cam], error = max((aa.x1+bb) / (cc.x1+dd)), cam starts from 1 */
    int i, j;
    double res_cam_aux[nCam][4];
    double res_cam[nCam][2];
    double List_cam[nCam];
    double List2[M];
    double max;
    
    for(i=0;i<M;i++)
    {
        a1[0] = a[i+0*M],  a1[1] = a[i+1*M],  a1[2] = a[i+2*M];
        c1[0] = c[i+0*M],  c1[1] = c[i+1*M],  c1[2] = c[i+2*M];
        b1 = b[i];
        d1 = d[i];
        reslist[i][0] = ((a1[0]*X[0] + a1[1]*X[1] + a1[2]*X[2]) + b1) / ((c1[0]*X[0] + c1[1]*X[1] + c1[2]*X[2]) + d1);
        reslist[i][1] = i+1;
        res_cam_aux[i%nCam][(int)(i/nCam)] = reslist[i][0];
    }
   
    for(i=0;i<nCam;i++)
    {
        max = 0;
        for(j=0;j<4;j++)
        {
            if( res_cam_aux[i][j] > max )
                max = res_cam_aux[i][j];
        }
        res_cam[i][0] = max;
        res_cam[i][1] = i+1;
    }

     
    qsort(res_cam,nCam,2*sizeof(double),compare_descent);
    
    *anchor = (int)res_cam[K-1][1];
    
    for(i=0;i<nCam;i++)
    {
        aid[(int)res_cam[i][1] - 1] = i+1;
        List[i] = (int)res_cam[i][1];
        
    }
}


double alpha_bf(double Ps[][3], int i)
{
    if(i<=0)
        return Ps[i][0]/2;
    else
        return (Ps[i][0] + Ps[i-1][0])/2;
}

double alpha_af(double Ps[][3], int i, int Ps_i)
{
    if(i>=Ps_i-1)
        return 1.01 * Ps[i][0];
    else
        return (Ps[i][0] + Ps[i+1][0])/2;
}


int tangent(double* a, double* b, double* c, double* d, \
            int M, double* X, double* dir, \
            double a_bf, double a_af, int cam1, int cam2)
{
    double err_1[2], err_2[2];
    double alphas[2] = {a_bf, a_af};
    get_error(a, b, c, d, X, dir, M, alphas, err_1, 2, cam1);
    get_error(a, b, c, d, X, dir, M, alphas, err_2, 2, cam2);

    if(  (err_1[0]>err_2[0] && err_1[1]>err_2[1]) || (err_1[0]<err_2[0] && err_1[1]<err_2[1])  )    
        return 1;
    else
        return 0;
}




double err_of_root(double* a1, double b1, double* c1, double d1, double* X, double* dir, double rt)
{
    double X1[3];
    X1[0] = X[0] + rt*dir[0]; 
    X1[1] = X[1] + rt*dir[1]; 
    X1[2] = X[2] + rt*dir[2];
    double v = ((a1[0]*X1[0] + a1[1]*X1[1] + a1[2]*X1[2]) + b1) / ((c1[0]*X1[0] + c1[1]*X1[1] + c1[2]*X1[2]) + d1);
    return v;
}


/* get error of camera cam for given list of alphas */
void get_error(double* a, double* b, double* c, double* d, \
               double* X, double* dir, \
               int M, double* alphas, double* errors, int n_alf, int cam)
{
    int nCam = ceil(M/4);
    int cid[4] = {cam-1, cam+nCam-1, cam+2*nCam-1, cam+3*nCam-1};
    double a0[3], a1[3], a2[3], a3[3]; 
    double b0, b1, b2, b3;    
    double c0[3], c1[3], c2[3], c3[3]; 
    double d0, d1, d2, d3;  
    double Xa[3];
    double r[4];
    int i,j;
    
    a0[0] = a[cid[0]+0*M],  a0[1] = a[cid[0]+1*M],  a0[2] = a[cid[0]+2*M];
    a1[0] = a[cid[1]+0*M],  a1[1] = a[cid[1]+1*M],  a1[2] = a[cid[1]+2*M];
    a2[0] = a[cid[2]+0*M],  a2[1] = a[cid[2]+1*M],  a2[2] = a[cid[2]+2*M];
    a3[0] = a[cid[3]+0*M],  a3[1] = a[cid[3]+1*M],  a3[2] = a[cid[3]+2*M];

    c0[0] = c[cid[0]+0*M],  c0[1] = c[cid[0]+1*M],  c0[2] = c[cid[0]+2*M];
    c1[0] = c[cid[1]+0*M],  c1[1] = c[cid[1]+1*M],  c1[2] = c[cid[1]+2*M];
    c2[0] = c[cid[2]+0*M],  c2[1] = c[cid[2]+1*M],  c2[2] = c[cid[2]+2*M];
    c3[0] = c[cid[3]+0*M],  c3[1] = c[cid[3]+1*M],  c3[2] = c[cid[3]+2*M];

    b0 = b[cid[0]];
    b1 = b[cid[1]];
    b2 = b[cid[2]];
    b3 = b[cid[3]];

    d0 = d[cid[0]];
    d1 = d[cid[1]];
    d2 = d[cid[2]];
    d3 = d[cid[3]];
    
    
    for(i=0;i<n_alf;i++)
    {
        Xa[0] = X[0] + alphas[i]*dir[0]; 
        Xa[1] = X[1] + alphas[i]*dir[1]; 
        Xa[2] = X[2] + alphas[i]*dir[2];
        
        r[0] = ((a0[0]*Xa[0] + a0[1]*Xa[1] + a0[2]*Xa[2]) + b0) / ((c0[0]*Xa[0] + c0[1]*Xa[1] + c0[2]*Xa[2]) + d0); 
        r[1] = ((a1[0]*Xa[0] + a1[1]*Xa[1] + a1[2]*Xa[2]) + b1) / ((c1[0]*Xa[0] + c1[1]*Xa[1] + c1[2]*Xa[2]) + d1); 
        r[2]= ((a2[0]*Xa[0] + a2[1]*Xa[1] + a2[2]*Xa[2]) + b2) / ((c2[0]*Xa[0] + c2[1]*Xa[1] + c2[2]*Xa[2]) + d2); 
        r[3] = ((a3[0]*Xa[0] + a3[1]*Xa[1] + a3[2]*Xa[2]) + b3) / ((c3[0]*Xa[0] + c3[1]*Xa[1] + c3[2]*Xa[2]) + d3); 
        
        double max = 0;
        for(j=0;j<4;j++)
        {
            if (r[j]>max)
                max = r[j];
        }
        errors[i] = max;
    }
}

/* get extremity */
void get_extreme(double* a, double* b, double* c, double* d, \
                   double lb, double ub, double* X, double* dir, \
                   int M, double Ps[][3], int* Ps_i)
{
    int cam,i,j;
    double alf_cands[12], values[12];
    int st;
    double roots[6], vs[6];
    int n_rt;
    double min;
    int nCam = ceil(M/4);

    for(cam=1;cam<=nCam;cam++) /* cam: cam id starting from 1 */
    {
        st = 0;
        int cid[4] = {cam-1, cam+nCam-1, cam+2*nCam-1, cam+3*nCam-1};
        for(i=0;i<3;i++)
        {
            int id_length = 4-i-1;
            n_rt = 0;
            f_root_positive( a, b, c, d, cid[i], (int*)(cid)+i+1, id_length, X, dir, lb, ub, M, roots, vs, &n_rt );
            
            for(j=0;j<n_rt;j++)
            {
                if(roots[j]>eps_root)
                {
                    alf_cands[st] = roots[j];
                    values[st] = vs[j];
                    st++;
                }
            }
        }
        min = Inf;
        int t = -1;
        for(i=0;i<st;i++)
        {
            if(values[i]<min)
            {
                min = values[i];
                t = i;
            }
        }
        if(t!=-1)
        {
            Ps[*Ps_i][0] = alf_cands[t];
            Ps[*Ps_i][1] = cam;
            Ps[*Ps_i][2] = 0;
            *Ps_i = *Ps_i + 1;
        }
    }
}
   
/* get intersections */
void get_intersection(double* a, double* b, double* c, double* d, \
                   double lb, double ub, double* X, double* dir, \
                   int M, double Ps[][3], int* Ps_i)
{
    double eps_root = 0.0000000001; /* 1e-10*/
    double eps_vid = 0.000001; /* 1e-6*/
    int n, m, i,j;
    double pair_sect[32][4];/*[alpha, value, cid1, cid2]*/
    double alphas[32], n_errors[32], m_errors[32];
    double roots[8], vs[8];
    int n_rt;
    int pid = 0; /* index of pair_sect */
    int nCam = ceil(M/4);
    
    for(n=0;n<nCam-1;n++)
    {
        
        int cid_1[4] = {n, n+nCam, n+2*nCam, n+3*nCam};
        for(m=n+1; m<nCam; m++)
        {
            int cid_2[4] = {m, m+nCam, m+2*nCam, m+3*nCam};
            pid = 0;
            for(i=0;i<4;i++)
            {
                f_root_positive( a, b, c, d, cid_1[i], cid_2, 4, X, dir, lb, ub, M, roots, vs, &n_rt );
                for(j=0;j<n_rt;j++)
                {
                    if(roots[j]>eps_root)
                    {
                        pair_sect[pid][0] = roots[j];
                        pair_sect[pid][1] = vs[j];
                        pair_sect[pid][2] = n+1;
                        pair_sect[pid][3] = m+1;
                        alphas[pid] = roots[j];
                        pid++;
                    }
                }
            }
            
            get_error( a, b, c, d, X, dir, M, alphas, n_errors, pid, n+1);
            get_error( a, b, c, d, X, dir, M, alphas, m_errors, pid, m+1);
            
            for(i=0;i<pid;i++)
            {
                if( n_errors[i]-pair_sect[i][1]<eps_vid && m_errors[i]-pair_sect[i][1]<eps_vid)
                {
                    Ps[*Ps_i][0] = pair_sect[i][0];
                    Ps[*Ps_i][1] = pair_sect[i][2];
                    Ps[*Ps_i][2] = pair_sect[i][3];
                    *Ps_i = *Ps_i + 1;
                }
            }
        }
    }
}

/* get positive roots of constraints */
void f_root_positive( double* a, double* b, double* c, double* d, \
                        int cid_i, int* cid_j, int id_length, \
                        double* X, double* dir, double lb, double ub, \
                        int M, double *roots, double* vs, int* n_rt )
{
    int i,j,J;
    double A,B,C;
    double temp;
    double rt;
    int pos;
    double a1[3], a2[3]; 
    double b1, b2;    
    double c1[3], c2[3]; 
    double d1, d2; 
    double v;
    
    a1[0] = a[cid_i+0*M],  a1[1] = a[cid_i+1*M],  a1[2] = a[cid_i+2*M];
    c1[0] = c[cid_i+0*M],  c1[1] = c[cid_i+1*M],  c1[2] = c[cid_i+2*M];
    b1= b[cid_i];
    d1= d[cid_i];

    *n_rt = 0;
    for(j=0;j<id_length;j++)
    {
        J = cid_j[j];
        a2[0] = a[J+0*M],  a2[1] = a[J+1*M],  a2[2] = a[J+2*M];
        c2[0] = c[J+0*M],  c2[1] = c[J+1*M],  c2[2] = c[J+2*M];
        b2= b[J];
        d2= d[J];

        A = (a1[0]*dir[0] + a1[1]*dir[1] + a1[2]*dir[2]) * (c2[0]*dir[0] + c2[1]*dir[1] + c2[2]*dir[2])  \
           -(a2[0]*dir[0] + a2[1]*dir[1] + a2[2]*dir[2]) * (c1[0]*dir[0] + c1[1]*dir[1] + c1[2]*dir[2]);

        B = (a1[0]*  X[0] + a1[1]*  X[1] + a1[2]*  X[2]) * (c2[0]*dir[0] + c2[1]*dir[1] + c2[2]*dir[2])  \
           +(a1[0]*dir[0] + a1[1]*dir[1] + a1[2]*dir[2]) * (c2[0]*  X[0] + c2[1]*  X[1] + c2[2]*  X[2])  \
           +(a1[0]*dir[0] + a1[1]*dir[1] + a1[2]*dir[2]) * d2                                            \
           +                                          b1 * (c2[0]*dir[0] + c2[1]*dir[1] + c2[2]*dir[2])  \
           -(a2[0]*  X[0] + a2[1]*  X[1] + a2[2]*  X[2]) * (c1[0]*dir[0] + c1[1]*dir[1] + c1[2]*dir[2])  \
           -(a2[0]*dir[0] + a2[1]*dir[1] + a2[2]*dir[2]) * (c1[0]*  X[0] + c1[1]*  X[1] + c1[2]*  X[2])  \
           -(a2[0]*dir[0] + a2[1]*dir[1] + a2[2]*dir[2]) * d1                                            \
           -                                          b2 * (c1[0]*dir[0] + c1[1]*dir[1] + c1[2]*dir[2]);

        C = (a1[0]*  X[0] + a1[1]*  X[1] + a1[2]*  X[2]) * (c2[0]*  X[0] + c2[1]*  X[1] + c2[2]*  X[2])  \
           +(a1[0]*  X[0] + a1[1]*  X[1] + a1[2]*  X[2]) * d2                                            \
           +                                          b1 * (c2[0]*  X[0] + c2[1]*  X[1] + c2[2]*  X[2])  \
           +                                          b1 * d2                                            \
           -(a2[0]*  X[0] + a2[1]*  X[1] + a2[2]*  X[2]) * (c1[0]*  X[0] + c1[1]*  X[1] + c1[2]*  X[2])  \
           -(a2[0]*  X[0] + a2[1]*  X[1] + a2[2]*  X[2]) * d1                                            \
           -                                          b2 * (c1[0]*  X[0] + c1[1]*  X[1] + c1[2]*  X[2])  \
           -                                          b2 * d1                                          ;

        temp = B*B-4*A*C;

        if(2*A==0) 
        {
            continue;
        }
        
        if(temp>0) /* 2 real roots*/
        {
            rt = (-B + sqrt(temp)) / (2*A);
            v = err_of_root(a1,b1,c1,d1,X,dir,rt);
            if(v>eps_root && rt>lb && rt>eps_root && rt<ub ) 
            {
                roots[*n_rt] = rt;
                vs[*n_rt] = v;
                *n_rt = *n_rt + 1;            
            }
            else
            {
                roots[*n_rt] = -1;
                vs[*n_rt] = -Inf;
                *n_rt = *n_rt + 1;            
            }
        
            rt = (-B - sqrt(temp)) / (2*A);
            v = err_of_root(a1,b1,c1,d1,X,dir,rt);

            if(v>eps_root && rt>lb && rt>eps_root && rt<ub ) 
            {
                roots[*n_rt] = rt;
                vs[*n_rt] = v;
                *n_rt = *n_rt + 1;                       
            }
            else
            {
                roots[*n_rt] = -1;
                vs[*n_rt] = -Inf;
                *n_rt = *n_rt + 1;            
            }
        }
        else if(temp==0) /* 1 real root */
        {
            rt = -B / (2*A);
            v = err_of_root(a1,b1,c1,d1,X,dir,rt);
            if(v>eps_root && rt>lb && rt>eps_root && rt<ub ) /*1e-15, valid root */
            {
                roots[*n_rt] = rt;
                vs[*n_rt] = v;
                *n_rt = *n_rt + 1;                       
            }
            else
            {
                roots[*n_rt] = -1;
                vs[*n_rt] = -Inf;
                *n_rt = *n_rt + 1;            
            }
        }
    }
}


