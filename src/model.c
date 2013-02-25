#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>

// #define DEBUG_MODEL

#ifdef DEBUG_MODEL
#define DEBUG(X) { X; }
#else
#define DEBUG(X)
#endif


/* Some debug print functions */
inline void PrintVector(const char *name, const double *v, const int *n)
{
    Rprintf("%s: [", name);
    if(*n>0)
    {
        Rprintf("%e", v[0]);
    }
    for(int c=1; c<*n; ++c)
    {
        Rprintf(" %e", v[c]);
    }
    Rprintf("]\n");
}

inline void PrintIntVector(const char *name, const int *v, const int *n)
{
    Rprintf("%s: [", name);
    if(*n>0)
    {
        Rprintf("%d", v[0]);
    }
    for(int c=1; c<*n; ++c)
    {
        Rprintf(" %d", v[c]);
    }
    Rprintf("]\n");
}

inline void PrintMatrix(const char *name, const double *v, const int *n)
{
    Rprintf("%s: [\n", name);
    for(int r=0; r<*n; ++r)
    {
        for(int c=0; c<*n; ++c)
        {
            Rprintf("%e ", v[c*(*n) + r]);
        }
        Rprintf("\n");
    }
    Rprintf("]\n");
}


inline void initNA(double *v, int n)
{
    /* Initialize the vector of doubles with NA */
    for(int x=0; x<n; x++) v[x] = NA_REAL;
}

void PrintParms(const int *n, 
                const double *K, 
                const double *a, 
                const double *q, 
                const double *d, 
                const double *W, 
                const int *producers, 
                const int *n_producers, 
                const int *consumers, 
                const int *n_consumers, 
                const double *rho, 
                const double *x, 
                const double *y, 
                const double *e, 
                const double *fe)
{
    /* Prints the parameters */
    Rprintf("n: [%d]\n", *n);
    Rprintf("K: [%e]\n", *K);
    PrintMatrix("a", a, n);
    Rprintf("q: [%e]\n", *q);
    PrintVector("d", d, n);
    PrintMatrix("W", W, n);
    PrintIntVector("producers", producers, n_producers);
    PrintIntVector("consumers", consumers, n_consumers);
    PrintVector("rho", rho, n);
    PrintVector("x", x, n);
    PrintMatrix("y", y, n);
    PrintMatrix("e", e, n);
    PrintMatrix("fe", fe, n);
}

void YodzisInnesState(const int *n_species,     /* n species */
                      const double *K,          /* Global carrying capacity */
                      const double *a,          /* n x n matrix */
                      const double *q,          /* Single number */
                      const double *d,          /* vector */
                      const double *W,          /* n x n matrix */
                      const int *producers,     /* vector */
                      const int *n_producers,   /* n items in producers */
                      const int *consumers,     /* vector */
                      const int *n_consumers,   /* n items in consumers */
                      const double *rho,        /* vector of length n */
                      const double *x,          /* vector of length n */
                      const double *y,          /* matrix of n x n */
                      const double *e,          /* n x n matrix */
                      const double *fe,         /* n x n matrix */
                      const double *B,          /* vector of length n */
                      double *dydt,             /* Output vector of length n */
                      double *growth,           /* Output vector of length n */
                      double *respiration,      /* Output vector of length n */
                      double *assimilation,     /* Output matrix of nxn */
                      double *consumption       /* Output matrix of nxn */
                      )
{
    /* Implements the model equations of Williams et al. 2007. Homage to 
       Yodzis and Innes 1992: scaling up feeding-based population dynamics to 
       complex ecological networks, in Rooney et al. (eds.), From Energetics 
       to Ecosystems: The Dynamics and Structure of Ecological Systems, 
       p 37–52. 
       Model equations are 2.17 and 2.18 on p 43. */

    /* Functional response and growth models are from Williams. 2008. Effects 
       of network and dynamical model structure on species persistence in 
       large model food webs. Theoretical Ecology. */

    const int n = *n_species;

    DEBUG(Rprintf("*************************\n"));
    DEBUG(PrintVector("B", B, n_species));
    DEBUG(PrintParms(n_species, K, a, q, d, W, producers, n_producers, 
                     consumers, n_consumers, rho, x, y, e, fe));

    /* Functional response numerator and denominator */
    /* Williams 2008 Eq 4 */
    double *fr_numerator = (double *) R_alloc(n*n, sizeof(double));
    initNA(fr_numerator, n*n);

    double *fr_denominator = (double *) R_alloc(n, sizeof(double));
    initNA(fr_denominator, n);
    if(TRUE)
    {
        const double exponent = 1.0 + *q;
        for(int consumer=0; consumer<*n_consumers; consumer++)
        {
            const int j = consumers[consumer];
            double col_sum = 0;
            for(int i=0; i<n; i++)
            {
                const int index = j*n + i;
                if(!ISNA(W[index]) && !ISNAN(W[index]))
                {
                    const double v = pow(B[i] / W[index], exponent);
                    fr_numerator[index] = v;
                    col_sum += v;
                }
            }

            fr_denominator[j] = 1.0 + d[j] * B[j] + col_sum;
        }
    }

    DEBUG(PrintMatrix("fr_numerator", fr_numerator, n_species));
    DEBUG(PrintVector("fr_denominator", fr_denominator, n_species));

    /* The assimilation matrix will hold the + term of equation 2.18 for each 
       spp */

    /* The consumption matrix will hold the - terms of equations 2.17 and 2.18 
       for each spp */

    /* Vectors that will hold the total gains from assimilation and 
       losses to consumption for each species.
       The assimilation_t vector is the +x[i]sum() term of equation 2.18 */
    double *assimilation_t = (double *) R_alloc(n, sizeof(double));
    memset(assimilation_t, 0, sizeof(double)*n);

    /* The consumption vector is the -sum() terms of equations 2.17 and 2.18 */
    double *consumption_t = (double *) R_alloc(n, sizeof(double));
    memset(consumption_t, 0, sizeof(double)*n);

    for(int col=0; col<n; col++)
    {
        for(int row=0; row<n; row++)
        {
            const int index = col*n + row;
            if(!ISNA(fr_numerator[index]) && !ISNAN(fr_numerator[index]) &&
               !ISNA(fr_denominator[col]) && !ISNAN(fr_denominator[col]))
            {
                const double v = x[col] * y[index] * B[col] * 
                                 fr_numerator[index] / fr_denominator[col];

                /* Store the values for this species in the matrices */
                assimilation[index] = v;
                consumption[index] = v / fe[index] / e[index];

                /* Add to the totals for each species */
                assimilation_t[col] += v;
                consumption_t[row] += consumption[index];
            }
        }
    }

    /* Producers - equation 2.17 using growth model Williams 2008 */
    for(int producer=0; producer<*n_producers; producer++)
    {
        const int j = producers[producer];
        double sum = 0;
        for(int competitor=0; competitor<*n_producers; competitor++)
        {
            const int k = producers[competitor];
            const int index = j + n*k;
            sum += a[index] * B[k];
        }

        growth[j] = rho[j] * B[j] * (1.0 - sum / *K);
        dydt[j] = growth[j] - consumption_t[j];
    }

    /* Consumers - equation 2.18 */
    for(int consumer=0; consumer<*n_consumers; consumer++)
    {
        const int j = consumers[consumer];
        respiration[j] = -x[j] * B[j];
        dydt[j] = respiration[j] + assimilation_t[j] - consumption_t[j];
    }

    DEBUG(PrintVector("growth", growth, n_species));
    DEBUG(Rprintf("\n"));

    DEBUG(PrintVector("respiration", respiration, n_species));
    DEBUG(Rprintf("\n"));

    DEBUG(PrintMatrix("consumption", consumption, n_species));
    DEBUG(PrintVector("consumption_t (row sums)", consumption_t, n_species));
    DEBUG(Rprintf("\n"));

    DEBUG(PrintMatrix("assimilation", assimilation, n_species));
    DEBUG(PrintVector("assimilation_t (col sums)", assimilation_t, n_species));
    DEBUG(Rprintf("\n"));

    DEBUG(PrintVector("dydt", dydt, n_species));

    DEBUG(Rprintf("*************************\n"));
    DEBUG(Rprintf("\n"));
}

void YodzisInnesFast(const int *neq,    /* n equations */ 
                     const double *t,   /* independent variable - time */
                     const double *B,   /* length *neq - current state */
                     double *dydt,      /* length *neq */
                     double *yout,      /* messy - see below */
                     const int *ip      /* messy - see below */
                     )
{
    /* WARNING - for reference only. Despite the name, this function is a bit 
       slower than the R wrapper method.

       Intended to be called directly from lsoda. Delegates work to 
       YodzisInnesState(). I thought that passing this C function to lsoda 
       would be quicker than using the R wrapper function YodzisInnesDyDt() as 
       the intermediate R step would be avoided. However, this function is 
       actually slower than the R wrapper function. 

       Using YodzisInnesDyDt():
            [1] "Simulation time:"
               user  system elapsed
              0.988   0.000   0.988

       Using YodzisInnesFast():
            [1] "Simulation time:"
               user  system elapsed
              0.980   0.084   1.122

        It might be slower because when using this scheme the output 
        parameters need discarded by RunChunk.LSODASimulation()
    */


    /* If you want to try this function, add to BuildModelParams():

   dll.ipar=c(length(producers), length(consumers), 
              producers.c, consumers.c), 
   dll.rpar=c(K,a,q,d,W,rho,x,y,e,fe), 
   dll.nout=2*NumberOfSpecies(community) + 
            2*NumberOfSpecies(community)^2)


       Add to RunChunk.LSODASimulation() after lsoda call:

    # Get rid of output variables. The 1+ is because first column is time
    chunk <- chunk[,1:(1+length(current.state))]


       Create the simulation using these parameters
    simulation <- LSODASimulation(model="YodzisInnesFast", 
                                  dllname='model', 
                                  ipar=params[['dll.ipar']],
                                  rpar=params[['dll.rpar']],
                                  nout=params[['dll.nout']],
    */

#ifdef DEBUG_MODEL
    const int debug_print_one=1;  // Used only in debug
#endif

    DEBUG(Rprintf("*************************\n"));
    DEBUG(PrintVector("t", t, &debug_print_one));
    DEBUG(PrintIntVector("neq", neq, &debug_print_one));

    if(neq<0) error("Must have at least 1 population");
    const int n = *neq;

    /* TASK 1 Extract integer parameters from ip */
    /* Format of ip: */
    /* 0    the number of output values */
    /* 1    the length out yout */
    /* 2    the length out ip */

    /* We expect */
    /* 3     number of producers */
    /* 4     number of consumers */
    /* 5...(5 + n producers)      indices of producers */
    /* (5 + n producers)...(5 + n producers + n consumers) indices of 
                                                           consumers */
    DEBUG(PrintIntVector("ip", ip, ip+2));

    if(ip[2]<5) error("Must have at least 5 integer parameters");
    const int *n_producers = ip + 3;
    const int *n_consumers = ip + 4;

    DEBUG(PrintIntVector("n_producers", n_producers, &debug_print_one));
    DEBUG(PrintIntVector("n_consumers", n_consumers, &debug_print_one));
 
    if(ip[2]!=(5+*n_producers+*n_consumers))
    {
        error("Unexpected number of integer parameters");
    }

    const int *producers = ip + 6;
    const int *consumers = producers + *n_producers;

    /* Check length of yout */
    const int n_outputs = n +       /* growth - vector of length n */
                          n +       /* respiration - vector of length n */
                          n*n +     /* assimilation - matrix of nxn */
                          n*n;      /* growth - vector of length n */

    const int n_params = 1 +        /* K - single value */
                         n*n +      /* a - matrix of nxn*/
                         1 +        /* q - single value */
                         n +        /* d - vector of length n */
                         n*n +      /* W - matrix of nxn*/
                         n +        /* rho - vector of length n */
                         n +        /* x - vector of length n */
                         n*n +      /* y - matrix of nxn*/
                         n*n +      /* e - matrix of nxn*/
                         n*n;       /* fe - matrix of nxn*/

    DEBUG(PrintIntVector("n_outputs", &n_outputs, &debug_print_one));
    DEBUG(PrintIntVector("n_params", &n_params, &debug_print_one));
    DEBUG(Rprintf("*************************\n"));

    if(ip[0]!=n_outputs)
    {
        error("Unexpected number of output parameters");
    }

    if(ip[1]!=n_outputs+n_params)
    {
        error("Unexpected yout length");
    }

    DEBUG(PrintVector("yout", yout, ip+1));
    DEBUG(Rprintf("*************************\n"));

    /* TASK 2 Extract double parameters from yout */
    double *growth=yout;                    /* Output vector of length n */
    double *respiration=growth+n;           /* Output vector of length n */
    double *assimilation=respiration+n;     /* Output matrix of nxn */
    double *consumption=assimilation+n*n;   /* Output matrix of nxn */

    const double *K=consumption+n*n;    /* Global carrying capacity */
    const double *a=K+1;                /* n x n matrix */
    const double *q=a+n*n;              /* Single number */
    const double *d=q+1;                /* vector of length n */
    const double *W=d+n;                /* n x n matrix */
    const double *rho=W+n*n;            /* vector of length n */
    const double *x=rho+n;              /* vector of length n */
    const double *y=x+n;                /* matrix of n x n */
    const double *e=y+n*n;              /* n x n matrix */
    const double *fe=e+n*n;             /* n x n matrix */

    YodzisInnesState(neq, 
                     K, 
                     a, 
                     q, 
                     d, 
                     W, 
                     producers, 
                     n_producers, 
                     consumers, 
                     n_consumers, 
                     rho, 
                     x, 
                     y, 
                     e, 
                     fe, 
                     B, 
                     dydt, 
                     growth, 
                     respiration, 
                     assimilation, 
                     consumption);
}
