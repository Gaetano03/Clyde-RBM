#include "write_Outputs.hpp"


void Config_stream ( prob_settings settings )
{

std::cout << std::endl;



}



void write_modes_sPOD ( const Eigen::MatrixXd &Phi_cut, 
                        const Eigen::MatrixXd &Coords, 
                        std::string flag_prob )
                        {

                            int Nr = Coords.rows(); 
                            std::string filename = "Modes_sPOD.dat";
                            std::ofstream flow_data;
                            flow_data.open(filename.c_str());

                            if ( flag_prob == "SCALAR" )
                            {

                                // Write row of Headers
                                flow_data << "\"PointID\"" << " ";
                                flow_data << "\"x\"" << " ";
                                flow_data << "\"y\"" << " ";

                                if ( Coords.cols() == 3 )
                                    flow_data << "\"z\"" << " ";

                                std::string phi;

                                for ( int i = 0; i < Phi_cut.cols(); i++ )
                                {

                                    phi = "\"Phi_" + std::to_string(i+1) + "\""; 
                                    flow_data << phi << " ";

                                }

                                flow_data << std::endl;

                                //Write fields
                                for ( int i = 0; i < Nr; i++ )
                                {

                                    flow_data << i+1 << " ";
                                    flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                                    if ( Coords.cols() == 3 )
                                        flow_data << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";
                                        

                                    for (int j = 0; j < Phi_cut.cols(); j++)
                                        flow_data << std::setprecision(12) << std::scientific << Phi_cut(i,j) <<  " ";           

                                flow_data << std::endl;

                                }
                            // Close file
                            flow_data.close();

                            } else if ( flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D" )
                            {

                                // Write row of Headers
                                flow_data << "\"PointID\"" << " ";
                                flow_data << "\"x\"" << " ";
                                flow_data << "\"y\"" << " ";

                                std::string phix;
                                std::string phiy;

                                for ( int i = 0; i < Phi_cut.cols(); i++ )
                                {

                                    phix = "\"Phi_x_" + std::to_string(i+1) + "\""; 
                                    flow_data << phix << " ";
                                    phiy = "\"Phi_y_" + std::to_string(i+1) + "\""; 
                                    flow_data << phiy << " ";

                                }

                                flow_data << std::endl;

                                //Write fields
                                for ( int i = 0; i < Nr; i++ )
                                {

                                    flow_data << i+1 << " ";
                                    flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";

                                    for (int j = 0; j < Phi_cut.cols(); j++)
                                    {

                                        flow_data << std::setprecision(12) << std::scientific << Phi_cut(i,j) <<  " ";
                                        flow_data << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j) <<  " ";            

                                    }

                                flow_data << std::endl;

                                }
                            // Close file
                            flow_data.close();

                            } else if ( flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D" )
                            {

                                // Write row of Headers
                                flow_data << "\"PointID\"" << " ";
                                flow_data << "\"x\"" << " ";
                                flow_data << "\"y\"" << " ";
                                flow_data << "\"z\"" << " ";

                                std::string phix;
                                std::string phiy;
                                std::string phiz;

                                for ( int i = 0; i < Phi_cut.cols(); i++ )
                                {

                                    phix = "\"Phi_x_" + std::to_string(i+1) + "\""; 
                                    flow_data << phix << " ";
                                    phiy = "\"Phi_y_" + std::to_string(i+1) + "\""; 
                                    flow_data << phiy << " ";
                                    phiy = "\"Phi_z_" + std::to_string(i+1) + "\""; 
                                    flow_data << phiz << " ";

                                }

                                flow_data << std::endl;

                                //Write fields
                                for ( int i = 0; i < Nr; i++ )
                                {

                                    flow_data << i+1 << " ";
                                    flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                                    flow_data << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";

                                    for (int j = 0; j < Phi_cut.cols(); j++)
                                    {

                                        flow_data << std::setprecision(12) << std::scientific << Phi_cut(i,j) <<  " ";
                                        flow_data << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j) <<  " ";
                                        flow_data << std::setprecision(12) << std::scientific << Phi_cut(2*Nr+i,j) <<  " ";            

                                    }

                                flow_data << std::endl;

                                }
                            // Close file
                            flow_data.close();


                            } else 
                            {

                                std::cout << "Set well problem flag! Exiting ..." << std::endl;
                                exit (EXIT_FAILURE);

                            }


                        }


void write_modes_DMD ( const Eigen::MatrixXcd &Phi_cut,
                    const Eigen::MatrixXd &Coords, 
                    std::string flag_prob )
                    {

                            int Nr = Coords.rows(); 
                            std::string filenameI = "Modes_DMD_Imag.dat";
                            std::string filenameR = "Modes_DMD_Real.dat";
                            std::ofstream flow_dataI, flow_dataR;

                            flow_dataI.open(filenameI.c_str());
                            flow_dataR.open(filenameR.c_str());

                            if ( flag_prob == "SCALAR" )
                            {

                                // Write row of Headers
                                flow_dataI << "\"PointID\"" << " ";
                                flow_dataI << "\"x\"" << " ";
                                flow_dataI << "\"y\"" << " ";

                                flow_dataR << "\"PointID\"" << " ";
                                flow_dataR << "\"x\"" << " ";
                                flow_dataR << "\"y\"" << " ";

                                if ( Coords.cols() == 3 )
                                {
                                    flow_dataI << "\"z\"" << " ";
                                    flow_dataR << "\"z\"" << " ";
                                }

                                std::string phi;

                                for ( int i = 0; i < Phi_cut.cols(); i++ )
                                {

                                    phi = "\"PhiI_" + std::to_string(i+1) + "\""; 
                                    flow_dataI << phi << " ";
                                    phi = "\"PhiR_" + std::to_string(i+1) + "\""; 
                                    flow_dataR << phi << " ";
                                }

                                flow_dataI << std::endl;
                                flow_dataR << std::endl;

                                //Write fields
                                for ( int i = 0; i < Nr; i++ )
                                {

                                    flow_dataI << i+1 << " ";
                                    flow_dataI << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_dataI << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                                    
                                    flow_dataR << i+1 << " ";
                                    flow_dataR << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_dataR << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";

                                    if ( Coords.cols() == 3 )
                                    {
                                        flow_dataI << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";
                                        flow_dataR << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";
                                    }    

                                    for (int j = 0; j < Phi_cut.cols(); j++)
                                    {
                                        flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(i,j).imag() <<  " ";
                                        flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(i,j).real() <<  " ";           
                                    }

                                flow_dataI << std::endl;
                                flow_dataR << std::endl;

                                }
                            // Close file
                            flow_dataI.close();
                            flow_dataR.close();

                            } else if ( flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D" )
                            {

                                // Write row of Headers
                                flow_dataI << "\"PointID\"" << " ";
                                flow_dataI << "\"x\"" << " ";
                                flow_dataI << "\"y\"" << " ";

                                flow_dataR << "\"PointID\"" << " ";
                                flow_dataR << "\"x\"" << " ";
                                flow_dataR << "\"y\"" << " ";

                                std::string phix;
                                std::string phiy;

                                for ( int i = 0; i < Phi_cut.cols(); i++ )
                                {

                                    phix = "\"PhiI_x_" + std::to_string(i+1) + "\""; 
                                    flow_dataI << phix << " ";
                                    phiy = "\"PhiI_y_" + std::to_string(i+1) + "\""; 
                                    flow_dataI << phiy << " ";
                                    phix = "\"PhiR_x_" + std::to_string(i+1) + "\""; 
                                    flow_dataR << phix << " ";
                                    phiy = "\"PhiR_y_" + std::to_string(i+1) + "\""; 
                                    flow_dataR << phiy << " ";

                                }

                                flow_dataI << std::endl;
                                flow_dataR << std::endl;

                                //Write fields
                                for ( int i = 0; i < Nr; i++ )
                                {

                                    flow_dataI << i+1 << " ";
                                    flow_dataI << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_dataI << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                                    flow_dataR << i+1 << " ";
                                    flow_dataR << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_dataR << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";

                                    for (int j = 0; j < Phi_cut.cols(); j++)
                                    {

                                        flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(i,j).imag() <<  " ";
                                        flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j).imag() <<  " ";
                                        flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(i,j).real() <<  " ";
                                        flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j).real() <<  " ";            

                                    }

                                flow_dataI << std::endl;
                                flow_dataR << std::endl;

                                }
                            // Close file
                            flow_dataI.close();
                            flow_dataR.close();

                            } else if ( flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D" )
                            {

                                // Write row of Headers
                                flow_dataI << "\"PointID\"" << " ";
                                flow_dataI << "\"x\"" << " ";
                                flow_dataI << "\"y\"" << " ";
                                flow_dataI << "\"z\"" << " ";

                                flow_dataR << "\"PointID\"" << " ";
                                flow_dataR << "\"x\"" << " ";
                                flow_dataR << "\"y\"" << " ";
                                flow_dataR << "\"z\"" << " ";

                                std::string phix;
                                std::string phiy;
                                std::string phiz;

                                for ( int i = 0; i < Phi_cut.cols(); i++ )
                                {

                                    phix = "\"PhiI_x_" + std::to_string(i+1) + "\""; 
                                    flow_dataI << phix << " ";
                                    phiy = "\"PhiI_y_" + std::to_string(i+1) + "\""; 
                                    flow_dataI << phiy << " ";
                                    phiy = "\"PhiI_z_" + std::to_string(i+1) + "\""; 
                                    flow_dataI << phiz << " ";
                                    phix = "\"PhiR_x_" + std::to_string(i+1) + "\""; 
                                    flow_dataR << phix << " ";
                                    phiy = "\"PhiR_y_" + std::to_string(i+1) + "\""; 
                                    flow_dataR << phiy << " ";
                                    phiy = "\"PhiR_z_" + std::to_string(i+1) + "\""; 
                                    flow_dataR << phiz << " ";

                                }

                                flow_dataI << std::endl;
                                flow_dataR << std::endl;

                                //Write fields
                                for ( int i = 0; i < Nr; i++ )
                                {

                                    flow_dataI << i+1 << " ";
                                    flow_dataI << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_dataI << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                                    flow_dataI << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";

                                    flow_dataR << i+1 << " ";
                                    flow_dataR << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_dataR << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                                    flow_dataR << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";

                                    for (int j = 0; j < Phi_cut.cols(); j++)
                                    {

                                        flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(i,j).imag() <<  " ";
                                        flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j).imag() <<  " ";
                                        flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(2*Nr+i,j).imag() <<  " ";

                                        flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(i,j).real() <<  " ";
                                        flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j).real() <<  " ";
                                        flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(2*Nr+i,j).real() <<  " ";            

                                    }

                                flow_dataI << std::endl;
                                flow_dataD << std::endl;

                                }
                            // Close file
                            flow_dataI.close();
                            flow_dataR.close();


                            } else 
                            {

                                std::cout << "Set well problem flag! Exiting ..." << std::endl;
                                exit (EXIT_FAILURE);

                            }

                    }



void write_coeffs_sPOD ( const Eigen::MatrixXd &Coeffs,
                        const std::vector<double> &t_vec,
                        const Eigen::VectorXd &lam )
                        {

                            int Ns = t_vec.size(); 

                            std::string filename = "Coeffs_sPOD.dat";
                            std::ofstream flow_data;
                            flow_data.open(filename.c_str());

                            // Write row of Headers
                            flow_data << "\"Time(s)\"" << " ";

                            std::string coef;

                            for ( int i = 0; i < Coeffs.cols(); i++ )
                            {

                                coef = "\"Coef_mode_" + std::to_string(i+1) + "\""; 
                                flow_data << coef << " ";

                            }                     

                            flow_data << "\"Lambda\"";
                            flow_data << std::endl;

                            // Write coefficients
                            for ( int i = 0; i < Coeffs.rows(); i++)
                            {

                                flow_data << std::setprecision(8) << t_vec[i] << " ";

                                for ( int j = 0; j < Coeffs.cols(); j++ )
                                    flow_data << std::setprecision(8) << Coeffs(i,j) << " ";

                                flow_data << std::setprecision(8) << lam(i);
                                flow_data << std::endl;

                            }

                            flow_data.close();


                        }


void write_TimeDynamics_DMD ( const Eigen::VectorXcd omega,
                            const Eigen::VectorXcd alfa,
                            const Eigen::VectorXd t)
                            {

                                std::ofstream flow_data;
                                flow_data.open(filename);

                                int Nt = t.size();
                                int Nm = alfa.size();

                                for ( int i = 0; i < Nt ; i ++)
                                    flow_data << std::setprecision(12) << std::scientific << t << " ";

                                flow_data << std::endl;

                                for ( int i = 0; i < Nm; i++ )
                                {
                                    for ( int j = 0; j < Nt; j++ )
                                        flow_data << std::setprecision(12) << std::scientific << alfa(i)*std::exp(omega(i)*t(j)) << " ";

                                    flow_data << std::endl;

                                }

                            }


void write_Reconstructed_fields ( Eigen::MatrixXd Rec,
                                    const Eigen::MatrixXd &Coords,
                                    std::string filename,
                                    std::string flag_prob )
                        {

                            int Nr = Coords.rows(); 

                            std::ofstream flow_data;
                            flow_data.open(filename);

                            if ( flag_prob == "SCALAR" )
                            {

                                // Write row of Headers
                                flow_data << "\"PointID\"" << " ";
                                flow_data << "\"x\"" << " ";
                                flow_data << "\"y\"" << " ";

                                if ( Coords.cols() == 3 )
                                    flow_data << "z" << " ";

 
                                flow_data << "\"Rec_field\"" << " ";
                                flow_data << std::endl;

                                //Write fields
                                for ( int i = 0; i < Nr; i++ )
                                {

                                    flow_data << i+1 << " ";
                                    flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                                    if ( Coords.cols() == 3 )
                                        flow_data << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";
                                        
                                    flow_data << std::setprecision(12) << std::scientific << Rec(i,0) <<  " ";           

                                flow_data << std::endl;

                                }
                            // Close file
                            flow_data.close();

                            } else if ( flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D" )
                            {

                                // Write row of Headers
                                flow_data << "\"PointID\"" << " ";
                                flow_data << "\"x\"" << " ";
                                flow_data << "\"y\"" << " ";

                                std::string Recx;
                                std::string Recy;

                                Recx = "\"Rec_x\""; 
                                flow_data << Recx << " ";
                                Recy = "\"Rec_y\""; 
                                flow_data << Recy << " ";
                                flow_data << std::endl;

                                //Write fields
                                for ( int i = 0; i < Nr; i++ )
                                {
                                    
                                    if ( abs(Rec(i,0)) < 1e-12 )    /*It*/
                                        Rec(i,0) = 0.0;             /*doesn't*/
                                    if ( abs(Rec(i,1)) < 1e-12 )    /*mean*/
                                        Rec(i,1) = 0.0;             /*anything*/
                                    flow_data << i+1 << " ";
                                    flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                                    flow_data << std::setprecision(12) << std::scientific << Rec(i,0) <<  " ";
                                    flow_data << std::setprecision(12) << std::scientific << Rec(i,1) <<  " ";            

                                    flow_data << std::endl;

                                }
                            // Close file
                            flow_data.close();

                            } else if ( flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D" )
                            {

                                // Write row of Headers
                                flow_data << "\"PointID\"" << " ";
                                flow_data << "\"x\"" << " ";
                                flow_data << "\"y\"" << " ";
                                flow_data << "\"z\"" << " ";

                                std::string Recx;
                                std::string Recy;
                                std::string Recz;

                                Recx = "\"Rec_x\""; 
                                flow_data << Recx << " ";
                                Recy = "\"Rec_y\""; 
                                flow_data << Recy << " ";
                                Recy = "\"Rec_z\""; 
                                flow_data << Recz << " ";
                                flow_data << std::endl;

                                //Write fields
                                for ( int i = 0; i < Nr; i++ )
                                {

                                    flow_data << i+1 << " ";
                                    flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                                    flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                                    flow_data << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";
                                    flow_data << std::setprecision(12) << std::scientific << Rec(i,0) <<  " ";
                                    flow_data << std::setprecision(12) << std::scientific << Rec(i,1) <<  " ";
                                    flow_data << std::setprecision(12) << std::scientific << Rec(i,2) <<  " ";            

                                flow_data << std::endl;

                                }
                            // Close file
                            flow_data.close();


                            } else 
                            {

                                std::cout << "Set well problem flag! Exiting ..." << std::endl;
                                exit (EXIT_FAILURE);

                            }


                        }
