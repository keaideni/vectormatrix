#ifndef CALCU_H
#define CALCU_H
#include "DMRG.h"
MatrixXd Kron(const MatrixXd& a, const MatrixXd& b);
void ReadTruncM(MatrixXd& A, const int& logo);

double ParticleNo(const Parameter& para);
double ParticleNo(const Parameter& para)
{
        Sub Sys, m(para,1);
        vector<MatrixXd> ParticleMat;
        ParticleMat.push_back(m.SysAdag()*m.SysA());

        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        ParticleMat.at(j)=(Kron(ParticleMat.at(j), m.SysEye()));
                        if(i==para.LatticeSize()/2)continue;
                        MatrixXd tempTrunc;
                        ReadTruncM(tempTrunc, i);
                        ParticleMat.at(j)=tempTrunc.adjoint()*ParticleMat.at(j)*tempTrunc;

                }
                Sys.Read(i-1);
                ParticleMat.push_back(Kron(Sys.SysEye(), m.SysAdag()*m.SysA()));
                if(i==para.LatticeSize()/2)break;
                MatrixXd tempTrunc;
                ReadTruncM(tempTrunc, i);
                ParticleMat.at(i-1)=tempTrunc.adjoint()*ParticleMat.at(i-1)*tempTrunc;
        }
        ParticleMat.push_back(m.SysAdag()*m.SysA());
        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        ParticleMat.at(para.LatticeSize()/2+j)=
                        (Kron(ParticleMat.at(para.LatticeSize()/2+j), m.SysEye()));
                        MatrixXd tempTrunc;
                        if(i==para.LatticeSize()/2)continue;
                        ReadTruncM(tempTrunc, para.LatticeSize()-i+1);
                        ParticleMat.at(para.LatticeSize()/2+j)
                        =tempTrunc.adjoint()*ParticleMat.at(para.LatticeSize()/2+j)*tempTrunc;

                }
                Sys.Read(para.LatticeSize()+2-i);
                ParticleMat.push_back(Kron(Sys.SysEye(), m.SysAdag()*m.SysA()));
                if(i==para.LatticeSize()/2)break;
                MatrixXd tempTrunc;
                ReadTruncM(tempTrunc, para.LatticeSize()-i+1);
                ParticleMat.at(para.LatticeSize()/2+i-1)
                =tempTrunc.adjoint()*ParticleMat.at(para.LatticeSize()/2+i-1)*tempTrunc;
        }


        

        MatrixXd wave;
        ReadTruncM(wave, 10000);

        vector<double> Particle;
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                Particle.push_back((wave.adjoint()*ParticleMat.at(i)*wave).trace());
        }for(int i=para.LatticeSize()/2; i<para.LatticeSize(); ++i)
        {
                Particle.push_back((wave*(ParticleMat.at(i)).transpose()*wave.adjoint()).trace());
        }
        ofstream outfile("./result/ParticleNo");outfile.precision(20);
        for(int i=0; i<para.LatticeSize(); ++i)
        {
                outfile<<"i= "<<(i+1)<<" ,ParticleNo= "<<Particle.at(i)<<endl;
        }
        outfile.close();

        return Particle.at(para.LatticeSize()/2);


}


double OrderParameter(const Parameter& para);
double OrderParameter(const Parameter& para)
{
        Sub Sys, m(para,1);
        vector<MatrixXd> OrderMat;
        OrderMat.push_back(m.SysA());

        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        OrderMat.at(j)=(Kron(OrderMat.at(j), m.SysEye()));
                        if(i==para.LatticeSize()/2)continue;
                        MatrixXd tempTrunc;
                        ReadTruncM(tempTrunc, i);
                        OrderMat.at(j)=tempTrunc.adjoint()*OrderMat.at(j)*tempTrunc;

                }
                Sys.Read(i-1);
                OrderMat.push_back(Kron(Sys.SysEye(), m.SysA()));
                if(i==para.LatticeSize()/2)break;
                MatrixXd tempTrunc;
                ReadTruncM(tempTrunc, i);
                OrderMat.at(i-1)=tempTrunc.adjoint()*OrderMat.at(i-1)*tempTrunc;
        }
        OrderMat.push_back(m.SysA());
        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        OrderMat.at(para.LatticeSize()/2+j)=
                        (Kron(OrderMat.at(para.LatticeSize()/2+j), m.SysEye()));
                        MatrixXd tempTrunc;
                        if(i==para.LatticeSize()/2)continue;
                        ReadTruncM(tempTrunc, para.LatticeSize()-i+1);
                        OrderMat.at(para.LatticeSize()/2+j)
                        =tempTrunc.adjoint()*OrderMat.at(para.LatticeSize()/2+j)*tempTrunc;

                }
                Sys.Read(para.LatticeSize()+2-i);
                OrderMat.push_back(Kron(Sys.SysEye(), m.SysA()));
                if(i==para.LatticeSize()/2)break;
                MatrixXd tempTrunc;
                ReadTruncM(tempTrunc, para.LatticeSize()-i+1);
                OrderMat.at(para.LatticeSize()/2+i-1)
                =tempTrunc.adjoint()*OrderMat.at(para.LatticeSize()/2+i-1)*tempTrunc;
        }


        

        MatrixXd wave;
        ReadTruncM(wave, 10000);

        vector<double> OrderNo;
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                OrderNo.push_back((wave.adjoint()*OrderMat.at(i)*wave).trace());
        }for(int i=para.LatticeSize()/2; i<para.LatticeSize(); ++i)
        {
                OrderNo.push_back((wave*(OrderMat.at(i)).transpose()*wave.adjoint()).trace());
        }
        ofstream outfile("./result/Order");outfile.precision(20);
        for(int i=0; i<para.LatticeSize(); ++i)
        {
                outfile<<"i= "<<(i+1)<<" ,Order= "<<OrderNo.at(i)<<endl;
        }
        outfile.close();

        return OrderNo.at(para.LatticeSize()/2);
}

void Correlation(const Parameter& para);
void Correlation(const Parameter& para)
{
        Sub m(para, 0);
        vector<MatrixXd> CorrMat;
        MatrixXd Adag1(m.SysAdag());

        CorrMat.push_back(m.SysAdag()*m.SysA());
        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        CorrMat.at(j)=Kron(CorrMat.at(j), m.SysEye());
                        
                        if(i==para.LatticeSize()/2)continue;
                        
                        MatrixXd tempTrunc;
                        ReadTruncM(tempTrunc, i);
                        CorrMat.at(j)=tempTrunc.adjoint()*CorrMat.at(j)*tempTrunc;
                        
                }
                CorrMat.push_back(Kron(Adag1, m.SysA()));
                Adag1=Kron(Adag1, m.SysEye());
                if(i==para.LatticeSize()/2)continue;
                MatrixXd tempTrunc;
                ReadTruncM(tempTrunc, i);
                CorrMat.at(i-1)=tempTrunc.adjoint()*CorrMat.at(i-1)*tempTrunc;
                Adag1=tempTrunc.adjoint()*Adag1*tempTrunc;
        }

        vector<double> Corr;
        MatrixXd wave;
        ReadTruncM(wave, 10000);
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                Corr.push_back((wave.transpose()*CorrMat.at(i)*wave).trace());
        }
        Sub Sys;
        Sys.Read(para.LatticeSize()/2+2);
        MatrixXd A(Kron(Sys.SysEye(), m.SysA()));
        Corr.push_back((wave.adjoint()*Adag1*wave*A.transpose()).trace());

        ofstream outfile("./result/Correlation");
        for(int i=0; i<Corr.size(); ++i)
        {
                outfile<<"R= "<<i<<" ,Corr(R)= "<<Corr.at(i)<<endl;
        }
        outfile.close();
}



MatrixXd Kron(const MatrixXd& a, const MatrixXd& b)
{
        MatrixXd ab(MatrixXd::Zero(a.rows()*b.rows(),a.cols()*b.cols()));

        int sizer(b.rows()),sizec(b.cols());
        for(int i=0; i<a.rows(); ++i)
        {
                for(int j=0; j<a.cols(); ++j)
                {
                        int startr(i*b.rows()),startc(j*b.cols());

                        ab.block(startr,startc,sizer, sizec)=a(i,j)*b;
                }
        }

        return ab;
}




#endif // CALCU_H
