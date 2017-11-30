#include "DMRG.h"
#include "test.h"
#include "Calcu.h"


int Sub::nmax;

int main(void)
{
        Parameter para;
        Sub::nmax=para.nmax();
        DMRG haha(para);
        
        cout<<haha.FEnergy()<<endl<<haha.Entropy()<<endl<<ParticleNo(para)<<endl;
}