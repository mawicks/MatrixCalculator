static long int currentflops = 0;
static long int currentsqrts = 0;
static long int totalflops = 0;
static long int totalsqrts = 0;
static char printflag = 0;

void addflops (n)
int n;
{
  currentflops += n;
}
void addsqrts (n)
int n;
{ currentsqrts += n;
}

void onflops ()
{
  printflag = 1;
}
void offflops ()
{
  printflag = 0;
}

void prtflops ()
{
  if (printflag) {
    printf ("Elapsed Flops: %ld Total Flops: %ld", currentflops,
	    totalflops += currentflops);
    printf (" Elapsed Sqrts: %ld Total Sqrts: %ld\n", currentsqrts,
	    totalsqrts += currentsqrts);
  }
  currentflops = 0;
  currentsqrts = 0;
}

