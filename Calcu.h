#ifndef CALCU_H
#define CALCU_H
void ReadTruncM(MatrixXd& A, const int& logo);

double ParticleNo(const Parameter& para);
double ParticleNo(const Parameter& para)
{
        Sub Sys, m(para,1);
        Sys.Read(para.LatticeSize()/2-1);

        MatrixXd wave;
        ReadTruncM(wave, 10000);

        Sub NewSys(para, Sys, m, 0);

        MatrixXd OP(NewSys.SysAdag()*NewSys.SysA());

        return (wave.adjoint()*OP*wave).trace();


}


double Order(const Parameter& para);
double Order(const Parameter& para)
{
        Sub Sys, m(para, 1);
        Sys.Read(para.LatticeSize()/2-1);

        MatrixXd wave;
        ReadTruncM(wave, 10000);
        Sub NewSys(para, Sys, m, 0);

        return (wave.adjoint()*NewSys.SysA()*wave).trace();
}








#endif // CALCU_H
