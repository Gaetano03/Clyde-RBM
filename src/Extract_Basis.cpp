#include "Extract_Basis.hpp"


int Nmod ( double En, Eigen::VectorXd K_pc )
{

    double En_content = 0.0;
    int count = 0;

    while ( En_content < En && count < K_pc.size() )
    {
        En_content = K_pc(count);
        count++;
    }

    return count;

}



void eig_sort( Eigen::VectorXd &lam, Eigen::MatrixXd &eig_vec ) {

    unsigned int swap_count = 1;
    double temp;
    Eigen::VectorXd temp_vec(eig_vec.rows());

    while (swap_count > 0)
    {

        swap_count = 0;

        for(unsigned int index = 1; index < lam.size(); index++)
        {

            if ( lam(index) > lam(index-1) )
            {

                temp = lam(index-1);
                lam(index-1) = lam(index);
                lam(index) = temp;

                temp_vec = eig_vec.col(index-1);
                eig_vec.col(index-1) = eig_vec.col(index);
                eig_vec.col(index) = temp_vec;

                swap_count++;

            }
        }
    }
}

// Should work only with singular values ( Not eigenvalues problems )
int SVHT ( Eigen::VectorXd lam, int m, int n )
{

    double median;
    int n_sv = lam.size();

    if ( n_sv%2 != 0 )
        median = lam((n_sv-1)/2);
    else
        median = 0.5*(lam(n_sv/2) + lam(n_sv/2-1));

    double beta = (double)std::min(m,n)/(double)std::max(n,m);
    double omega = 0.56*std::pow(beta,3.0) - 0.95*std::pow(beta,2.0) 
                    + 1.82*beta + 1.43;
    double tau = omega*median;

    double eps = 1000;
    int Nm = 0;

    do
    {
        eps = lam(Nm);
        Nm++;
    } while ( eps > tau );

    return Nm-1;

}


Eigen::MatrixXd SPOD_basis( const Eigen::MatrixXd &snap_set,
                                Eigen::VectorXd &lam,
                                Eigen::VectorXd &K_pc,
                                Eigen::MatrixXd &eig_vec,
                                const int Nf,    
                                std::string bc_flag , 
                                std::string filter_flag ,  
                                double sigma )
{

    int count;
    int Nr = snap_set.rows();
    int Ns = snap_set.cols();

    Eigen::MatrixXd R_f(Ns, Ns);
    Eigen::MatrixXd phi_c(Nr, Ns);
    Eigen::VectorXd mean(Nr);
    Eigen::VectorXd g(2*Nf+1);

    R_f.setZero(Ns, Ns);

    if ( Nf == 0) 
    {                           //Calculating R POD
        R_f = snap_set.transpose()*snap_set;

    } else
    {                                      //Calculating R SPOD
    
        if ( filter_flag == "BOX" )
        {             //with BOX filter  

            for (int k = 0; k <= 2*Nf; k++ )
                g(k) = 1.0/(2.0*(double)Nf+1.0); 

        } else if ( filter_flag == "GAUSSIAN" )
        {             //with GAUSSIAN filter
            
            if ( sigma == 0.0)
            {
                std::cout << "sigma = 0 then only POD could be performed" << std::endl;
                exit (EXIT_FAILURE);
            }

            for ( int k = 0; k <= 2*Nf; k++ )
                g(k) = exp(-k*k/(2.0*sigma*sigma));
            

        } else
        {
            std::cout << "Filter selected not available" << std::endl;
            exit (EXIT_FAILURE);
        }

        Eigen::MatrixXd R = snap_set.transpose()*snap_set;

        //Zero-padded boundary conditions
        if (bc_flag == "ZERO"){

            for ( int i = 0; i < Ns; i++ )
            {
                for ( int j = 0; j < Ns; j++ )
                {

                    count = 0;

                    for ( int k = -Nf; k <= Nf; k++ ){
                        if ( i + k < 0 || j + k < 0 || i + k >= Ns || j + k >= Ns )
                            R_f(i,j) += 0.0;
                        else
                            R_f(i,j) += g(count)*R(i+k, j+k);
                        
                        count++; 
        
                    }
                }
            }
        } else 
        {
            std::cout << "Boundary condition not implemented " << std::endl;
            exit (EXIT_FAILURE);
        }
    }


    Eigen::EigenSolver<Eigen::MatrixXd> es(R_f); 
    lam = es.eigenvalues().real();
    eig_vec = es.eigenvectors().real();
    eig_sort( lam, eig_vec);

    double sum = 0;

    for (int i = 0; i < Ns; i++){
        sum += lam(i)/lam.sum();
        K_pc(i) = sum;
    }

    // double tol = lam(0)*1e-12;
    double tol = 1e-16;
    phi_c = snap_set*eig_vec;

    count = 0;
    while ( count < lam.size() && std::abs(lam(count)) > tol)
            count++;

    Eigen::MatrixXd phi(Nr,count);
    for ( int i = 0 ; i < count ; i++ )
        phi.col(i) = phi_c.col(i)/sqrt(lam(i));


    return phi;   


}



Eigen::MatrixXcd DMD_basis ( const Eigen::MatrixXd &snap_set,
                            Eigen::VectorXcd &lam,
                            Eigen::MatrixXcd &eig_vec,
                            Eigen::VectorXd &lam_POD,
                            Eigen::MatrixXd &eig_vec_POD,
                            const int r )
{   

    int Ns = snap_set.cols() - 1;
    int Nm;
    //                                Eigen::MatrixXd U = SPOD_basis(snap_set.leftCols(Ns), lam_POD, K_pc, eig_vec_POD );
    Eigen::BDCSVD<Eigen::MatrixXd> svd( snap_set.leftCols(Ns), 
                                        Eigen::ComputeThinU | Eigen::ComputeThinV );
    lam_POD = svd.singularValues();
    eig_vec_POD = svd.matrixV();
    eig_sort(lam_POD, eig_vec_POD);

    Eigen::MatrixXd U = svd.matrixU();                         
    Eigen::MatrixXd Sig_inv = Eigen::MatrixXd::Zero(U.cols(), U.cols());
    //std::cout << "Number of non-zero modes : " << U.cols() << std::endl;

    // for ( int i = 0; i < U.cols(); i++ )
        // Sig_inv(i, i) = 1.0/std::sqrt(lam_POD(i)); 

    for ( int i = 0; i < U.cols(); i++ )
        Sig_inv(i, i) = 1.0/lam_POD(i); 

    // int Nm = Nmod( En, K_pc );

    // if ( En == 1.0 )
    //     Nm = U.cols();

    if ( r == 0)
    {
        Nm = SVHT ( lam_POD, Ns, snap_set.rows() );
        std::cout << "DMD-rank from SVHT : " << Nm << std::endl;
    }
    else
    {                    
        Nm = std::min(r, Ns);
        std::cout << "DMD user-defined rank : " << Nm << std::endl;
    }

    Eigen::MatrixXd Atilde = U.leftCols(Nm).transpose()*snap_set.rightCols(Ns)*
                                eig_vec_POD.leftCols(Nm)*Sig_inv.block(0,0,Nm,Nm);


    //Eigen::VectorXcd Full_lam;

    if ( Atilde.size() == 1 && Atilde(0,0) == 0.0 )
    {
        lam = Eigen::VectorXcd::Zero(1);
        eig_vec = Eigen::MatrixXcd::Ones(1,1);
    }
    else 
    {
        Eigen::EigenSolver<Eigen::MatrixXd> es(Atilde); 
        // Full_lam = es.eigenvalues();
        lam = es.eigenvalues();
        eig_vec = es.eigenvectors();
    }

    Eigen::MatrixXcd appo = snap_set.rightCols(Ns)*
                            eig_vec_POD.leftCols(Nm)*Sig_inv.block(0,0,Nm,Nm);                              

    
    //This is not Working well!!!
    //-----------Return only modes with frequencies : f < f_sampling/10 with f_sampling reading snapshots frequency------
    //Attenzione!!! Anche in tal caso la matrice di autovettori  rimane sempre la stessa (N_rank_DMDxN_rank_DMD)
    //quindi il calcolo dei coefficienti va effettuato tramite la funzione Calculate_Coefficients_DMD_Exact
    // Eigen::MatrixXcd Phi = appo*eig_vec;
    // std::complex<double> omegaj;

    // std::vector<int> idx_ny_modes = {};
    // std::cout << "Full Lam : \n" << Full_lam << std::endl;
    // for ( int i = 0; i < Full_lam.size(); i++ )
    // {
    //     omegaj = std::log(Full_lam(i));

    //     if ( (std::abs(omegaj.imag())) <= 1.0/8.0 )
    //         idx_ny_modes.push_back(i);
    // }

    // Eigen::MatrixXcd Phi_nyq(snap_set.rows(),idx_ny_modes.size());
    // lam = Eigen::VectorXcd::Zero(idx_ny_modes.size());

    // for ( int i = 0; i < idx_ny_modes.size(); i++ )
    // {
    //     Phi_nyq.col(i) = Phi.col(idx_ny_modes[i]);
    //     lam(i) = Full_lam(idx_ny_modes[i]);
    // } 

    // std::cout << "Final number of non-spurious modes (f < f_sampling/10) : " << Phi_nyq.cols() << std::endl;

    // return Phi_nyq;


    //This is the only choice for now
    //------------Return Modes considering also spurious ones (f > Nyquist frequency)----------------
    // Non divido per lambda
    //Se non divido per lambda mi trovo un timestep avanti (perchè? non so ancora)
    return appo*eig_vec;

    // Divido per lambda

    //                                 Eigen::MatrixXcd phi(snap_set.rows(), Nm);
    //                                 for (int i = 0; i < Nm; i++)
    //                                     phi.col(i) = 1.0/lam(i)*appo*eig_vec.col(i);
    // // 
    //                                 return phi;

    //Standard DMD
    // return U*eig_vec;
}


//The vector of lam and eig_vec DMD has to contain only the selected modes for reconstruction
Eigen::VectorXcd Calculate_Coefs_DMD ( const Eigen::MatrixXcd &eig_vec,
                                    const Eigen::MatrixXcd &eig_vec_POD,
                                    const Eigen::VectorXcd &lam,
                                    const Eigen::VectorXcd &lam_POD,
                                    const int Ns)   //Number of original snapshots minus 1
{
    
    int Nm = lam.size();
    Eigen::MatrixXcd Sigma = Eigen::MatrixXcd::Zero(Nm, Nm);
    Eigen::MatrixXcd V_and(Nm, Ns); 

    for ( int i = 0; i < Nm; i++)
        Sigma(i, i) = lam_POD(i);

    for ( int i = 0; i < Nm; i++ )
    {
        for ( int j = 0; j < Ns; j++ )
            V_and(i,j) = std::pow(lam(i), j);                                                                                         
    }

    Eigen::MatrixXcd Y_sq = eig_vec.conjugate().transpose()*eig_vec;
    Eigen::MatrixXcd V_and_sq = V_and*V_and.conjugate().transpose();
    V_and_sq = V_and_sq.conjugate();

    Eigen::MatrixXcd M1(Nm, Nm);

    for ( int i = 0; i < Nm; i++ )
    {

        for ( int j = 0; j < Nm; j++)
            M1(i,j) = Y_sq(i,j)*V_and_sq(i,j);

    } 

    Eigen::MatrixXcd M2 = V_and*eig_vec_POD.leftCols(Nm)*Sigma.conjugate()*eig_vec;
    Eigen::VectorXcd dum = M2.diagonal();
    dum = dum.conjugate();

    return M1.inverse()*dum;
    

} 


Eigen::VectorXcd Calculate_Coefs_DMD_exact ( const Eigen::MatrixXd &sn_set,  //matrix of first Ns-1 snaps 
                                            const Eigen::VectorXcd &lam,  //slow eigenvalues
                                            const Eigen::MatrixXcd &Phi ) //slow exact DMD modes
{

    int Nm = lam.size();
    int Ns = sn_set.cols();
    Eigen::MatrixXcd V_and(Nm, Ns);

    for ( int i = 0; i < Nm; i++ )
    {
        for ( int j = 0; j < Ns; j++ )
            V_and(i,j) = std::pow(lam(i), j);                                                                                         
    }

    Eigen::MatrixXcd Phi_sq = Phi.transpose().conjugate()*Phi;
    Eigen::MatrixXcd V_sq = V_and*V_and.conjugate().transpose();
    Eigen::MatrixXcd V_sqc = V_sq.conjugate();

    Eigen::MatrixXcd P = Phi_sq.cwiseProduct(V_sqc);

    Eigen::MatrixXcd appo = V_and*sn_set.transpose()*Phi;
    Eigen::VectorXcd dum = appo.diagonal();
    Eigen::VectorXcd dumc = dum.conjugate();

    return P.inverse()*dumc;

}



std::vector<node_mrDMD> mrDMD_basis( Eigen::MatrixXd &snap_set,
                                    std::vector<node_mrDMD> &nodes,                          
                                    const int r,
                                    double dts,                                                
                                    double t_0,
                                    int level,
                                    int bin_num,
                                    int offset,
                                    int max_levels,
                                    int max_cycles)
{
    std::cout << "--------LEVEL " << level << "---------------" << std::endl << std::endl;                                        
    const double PI = 3.1415926535;
    int nyq = 8*max_cycles;                         //number of snaps needed to capture cycles

    int N = snap_set.rows();                        //dimension of each snapshot
    int bin_size = snap_set.cols();                 //number of total snapshots available for the particular time bin

    if ( bin_size < nyq )
    {
        std::cout << "Max resolution possible reached ..." <<
                    "\n Returning to above level" << std::endl << std::endl;
        return nodes;
    }
    int step = std::floor( bin_size/nyq );
    Eigen::MatrixXd _snap_set(N, nyq);

    for ( int i = 0, jj = 0; i < nyq; i += step, jj++ )
        _snap_set.col(jj) = snap_set.col(i);

    Eigen::VectorXd lam_POD;
    Eigen::VectorXcd lam_DMD;
    Eigen::MatrixXd eig_vec_POD;
    Eigen::MatrixXcd eig_vec_DMD;

    //Perform pure DMD in the time bin
    Eigen::MatrixXcd Phi = DMD_basis( _snap_set,
                                        lam_DMD,
                                        eig_vec_DMD,
                                        lam_POD,
                                        eig_vec_POD,
                                        r );

    std::cout << " Singular values : " << lam_POD << std::endl;
    //select only the slow modes (only if max_level > 0 (not pure DMD))
    double rho = (double)max_cycles/(double)bin_size;
    if ( max_levels == 0 )
        rho = 1e6;

    std::vector<int> slow_idx;

    for ( int i =0; i < lam_DMD.size(); i++ )
    {

        if ( (std::abs(std::log(lam_DMD(i))/(2*PI*step))) <= rho )  //define PI
            slow_idx.push_back(i);
    }

    int n = slow_idx.size();

    
    Eigen::MatrixXcd Modes(N,n);
    Eigen::VectorXcd lam_slw(n);
    Eigen::MatrixXcd eig_vec_slw(eig_vec_DMD.rows(), n);
    Eigen::VectorXcd b_opt;    
    Eigen::MatrixXcd Psi(n, bin_size);           //time evolution matrix
    Eigen::MatrixXcd D_dmd;

    //Calculate the correspondent coefficients and time evolution
    if ( n > 0 )
    {

        Eigen::MatrixXcd V_and(lam_slw.size(), bin_size);

        for ( int i = 0; i < n ; i++ )
        {

            Modes.col(i) = Phi.col(slow_idx[i]);                                                
            lam_slw(i) = lam_DMD(slow_idx[i]);
            eig_vec_slw.col(i) = eig_vec_DMD.col(slow_idx[i]);

        }

        b_opt = Calculate_Coefs_DMD_exact ( snap_set.leftCols(bin_size - 1),
                                            lam_slw,
                                            Modes );

        for ( int i = 0; i < lam_slw.size(); i++ )
        {
            for ( int j = 0; j < bin_size; j++ )
                V_and(i,j) = std::pow(std::pow(lam_slw(i), 1/step), j);                                                                                         
        }                                                                                    
//-----------Remove after check----------------------------------------------------------
// std::cout << " size V_and : [" << V_and.rows() << ", " << V_and.cols() << "]" << std::endl;
// std::cout << " size b_opt : [" << b_opt.rows() << ", " << b_opt.cols() << "]" << std::endl;

        for ( int i = 0; i < bin_size; i++ )
        {
            Psi.col(i) = b_opt.cwiseProduct(V_and.col(i));
        }
        D_dmd = Modes*Psi;
        
//------------Remove after check-------------------------------------------------------------
// std::cout << " size Modes : [" << Modes.rows() << ", " << Modes.cols() << "]" << std::endl;
// std::cout << " size Psi : [" << Psi.rows() << ", " << Psi.cols() << "]" << std::endl;                                     
// std::cout << " size D_dmd : [" << D_dmd.rows() << ", " << D_dmd.cols() << "]" << std::endl;
    }
    else        //initialize b_opt  and Psi as empty vector and matrix
    {

        b_opt = Eigen::VectorXcd::Zero(0);
        Psi = Eigen::MatrixXcd::Zero(0,0);
        D_dmd = Eigen::MatrixXcd::Zero(N,bin_size);

    }

    //Subtracting influence of slow modes
    snap_set = snap_set - D_dmd.real();

    //Storing all the necessary information in the node
    node_mrDMD node;
    node.l = level;
    node.bin_num = bin_num;
    node.bin_size = bin_size;      
    node.start = offset;           
    node.stop = offset + bin_size - 1; 
    node.step = step;          
    node.rho = rho;                
    node.r = lam_DMD.size();                    
    node.n = n;                    
    node.t_begin = t_0;
    node.t_end = t_0 + (bin_size - 1)*dts;
    node.dt = dts*step;
    node.lam = lam_slw;                  
    node.Modes = Modes;                
    node.Psi = Psi;         //time-dynamics transpose               
    node.Coefs = b_opt; 

    //Output nodes information
    std::cout << "----> Node  Info : " << std::endl << std::endl;
    std::cout <<  "Time interval : [" << node.t_begin << 
                    ", " << node.t_end << "]" << std::endl;    
    std::cout << " Snapshot interval (snaps index) : " << 
                    node.start << "\t" << node.stop << std::endl;
    std::cout << " bin_size : " << node.bin_size << std::endl;
    std::cout << " bin_num : " << node.bin_num << std::endl;
    std::cout << " DMD-rank : " << node.r << std::endl;
    std::cout << " Number of slow modes : " << node.n << std::endl;
    std::cout << std::endl << std::endl;

    nodes.push_back(node);

    //Apply DMD recursively
    int split = ceil(bin_size/2);

    if ( level < max_levels && split > 1 )
    {    

        Eigen::MatrixXd snap_set1;
        Eigen::MatrixXd snap_set2;
        
        if ( bin_size%2 != 0 )
        {
            snap_set1 = snap_set.leftCols(split + 1);
            nodes = mrDMD_basis( snap_set1,
                                    nodes,                         
                                    r,
                                    dts,
                                    t_0,
                                    level+1,
                                    2*bin_num,
                                    offset,
                                    max_levels,
                                    max_cycles);
            snap_set2 = snap_set.rightCols(split + 1);
            nodes = mrDMD_basis( snap_set2, 
                                nodes,                         
                                r,
                                dts,
                                t_0 + (double)split*dts,
                                level+1,
                                2*bin_num+1,
                                offset+split,
                                max_levels,
                                max_cycles);

        }
        else
        {
            
            snap_set1 = snap_set.leftCols(split + 1);
            nodes = mrDMD_basis( snap_set1,
                                    nodes,                         
                                    r,
                                    dts,
                                    t_0,
                                    level+1,
                                    2*bin_num,
                                    offset,
                                    max_levels,
                                    max_cycles);

            snap_set2 = snap_set.rightCols(split);
            nodes = mrDMD_basis( snap_set2,
                                nodes,                         
                                r,
                                dts,
                                t_0 + (double)(split)*dts,
                                level+1,
                                2*bin_num+1,
                                offset+split,
                                max_levels,
                                max_cycles);

        }
    }

    return nodes;

}





