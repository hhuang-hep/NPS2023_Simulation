#include "TSystem.h"
#include "TMatrix.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TCutG.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TStyle.h"
#include <iostream>
#include <sstream>
#include "TTreeIndex.h"
#include "TChainIndex.h"
#include "TCaloBase.h"
#include "TCaloEvent.h"
#include "TCaloGeometry.h"
#include "TDVCSDB.h"
#include "TDVCSEvent.h"
using namespace std;

void GetRunRange(const string &kinc_param, vector<int> &minrun, vector<int> &maxrun);
void GetRunNumbersAndCharge(int target_flag, vector<int> minrun, vector<int> maxrun, vector <int> &run_number, vector <double> &total_charge);
double GetVertexZ(double GV_x, double HMS_Angle, double hsytar, double hsyptar); // Reconstruction of vertex z
Int_t bnConv_OldToNew(int ibn_old); // Convert the PMT number to the current version
Int_t bnConv_NewToOld(int ibn_new); // Convert the PMT number to the simulation version
TDVCSDB *db = new TDVCSDB("dvcs", "clrlpc", 3306, "hhuang", "");

void Reconstruction(string kinc_param, int target_flag = 0, int i_job = 1, int gen_type = 0)
{   
    gSystem->Load("libDVCS.so");

    // Specify the target type
    Double_t M_targ = 0;
    TString Tar = "";
    if(target_flag == 0){ // proton
        M_targ = 0.938272;
        Tar = "LH2";
    }
    else if(target_flag == 1){ // neutron
        M_targ = 0.939565;
        Tar = "LD2";
    }
    else{
        cout << "Error: Invalid target flag!" << endl; 
        return;
    }

    // 2x2 clustering threshold in GeV
    Double_t clusThr;
    if(gen_type == 0) clusThr = 0.5; // for DVCS photons
    else if(gen_type == 1) clusThr = 0.2; // for pi0 photons
    else {
        cout << "Error: invalid reaction type!" << endl;
        return;
    }

    // The minimum and maxrun numbers of the kinematics
    vector<int> minrun, maxrun;
    GetRunRange(kinc_param, minrun, maxrun);
    // cout<< minrun.size()<<endl;
    // for(int i = 0; i<minrun.size(); i++) cout<<minrun.at(i)<<" "<<maxrun.at(i)<<endl;

    // Total charge of the production runs of this kinematics and target
    vector <int> run_number;
    vector <double> total_charge_frac;
    GetRunNumbersAndCharge(target_flag, minrun, maxrun, run_number, total_charge_frac);

    // Get the run information from the database
    Double_t HMS_angle = *db->GetEntry_d("SIMU_param_HMSangle", run_number[0]);
    Double_t caloDist = *db->GetEntry_d("CALO_geom_Dist", run_number[0]);
    Double_t caloAngle = *db->GetEntry_d("CALO_geom_Yaw", run_number[0]);
    Int_t *caloMaskBlock_temp = new Int_t[1080]; // Get the mask block information in NPS numbering Scheme
    Int_t *caloMaskBlock = new Int_t[1080]; // simulation numbering scheme
    caloMaskBlock_temp = db->GetEntry_i("CALO_flag_MaskBlock", run_number[0]);
    // Convert the mask block information to the simulation numbering scheme
    for(Int_t i = 0; i < 1080; i++){
        caloMaskBlock[i] = caloMaskBlock_temp[bnConv_OldToNew(i)];
    }
    // for(int i = 0; i < 1080; i++) cout << caloMaskBlock[i] << " ";
    // cout << endl;

    string gen_type_str = "";
    if(gen_type == 0) gen_type_str = "DVCS";
    else if(gen_type == 1) gen_type_str = "pi0";
    else {cout<<"Error: invalid reaction type!"<<endl; return;}
    cout<<"Start Processing photon reconstruction of "<< gen_type_str.c_str()<<" for "<<kinc_param.c_str()<<endl;
    cout<<"========== Kinematic information =========="<<endl;
    cout<<"Target type: ";
    if(target_flag == 0) cout<<"LH2"<< endl;
    if(target_flag == 1) cout<<"LD2"<< endl;
    cout<<"HMS angle: "<<HMS_angle*TMath::RadToDeg()<<" degrees"<<endl;
    cout<<"NPS distance: "<<caloDist<<" cm"<<endl;
    cout<<"NPS angle: "<<caloAngle*TMath::RadToDeg()<<" degrees"<<endl;
    cout<<"Number of production runs: "<<run_number.size()<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;
    cout<<"========== NPS mask blocks =========="<<endl;
    for(int i = 0; i < 1080; i++) if(caloMaskBlock_temp[i] == 1) cout<<i<<" "; // output in NPS numbering scheme
    cout<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;
    cout<<"========== Total charge to account for the dead blocks=========="<<endl;
    for(int irun = 0; irun < run_number.size(); irun++){
        cout<<"Run: "<<run_number[irun]<<", Total charge: "<<total_charge_frac[irun]<<" uC"<<endl;
    }
    cout<<endl;

    // Fraction of accumulated charge after adding each run
    Double_t sum_charge = 0;
    for(int irun = 0; irun < run_number.size(); irun++) sum_charge += total_charge_frac[irun];
    for(int irun = 0; irun < run_number.size(); irun++) total_charge_frac[irun] = total_charge_frac[irun]/sum_charge;
    for(int irun = 0; irun < run_number.size(); irun++) total_charge_frac[irun+1] += total_charge_frac[irun]; // cumulative charge fraction
    for(int irun = 0; irun < run_number.size(); irun++) cout<<total_charge_frac[irun]<<endl; // cumulative charge fraction

    // Input file from the Geant4 simulation
    TString g4_inputDir = "generated_events";
    TString g4_filename = Form("%s/%s_%s_%d.root", g4_inputDir.Data(), kinc_param.c_str(), Tar.Data(), i_job);
    // TString g4_filename = "nps_hms_simulation_x36_5_3_LH2_1.root";
    
    // Event generator and Geant4 NPS simulation______________________________________
    TFile *infile_g4 = TFile::Open(Form("%s", g4_filename.Data()));
    if(!infile_g4){
        cout<<"Error: can not open "<<g4_filename.Data()<<endl;
        return;
    }
    TTree *t_g4 = (TTree*)infile_g4->Get("t_dvcs");

    Int_t evtNb; // Event number from event generator
    Double_t GV_x,GV_y, GV_z; // Generated vertex position, the raster and beam offsets are accounted
    Double_t GQ2, GxB, Gt, Gphi, psf; // Q2, xB, t, phi and phase space factor from DVCS gen
    Double_t RIE_px, RIE_py, RIE_pz; // Beam energy of the kinematics
    Double_t GIE_px, GIE_py, GIE_pz; // Beam energy after the pre-vertex external radiation correction
    Double_t GSE_px, GSE_py, GSE_pz; // Generated scattered electron with the internal radiation correction
    Double_t RSE_px, RSE_py, RSE_pz; // Reconstructed scattered electon from HMS MC simulation with external radiation correction
    Double_t GP_px, GP_py, GP_pz; // Generated real photon from DVCS gen

    Double_t GP_px1, GP_py1, GP_pz1; // Generated real photon from pi0 gen
    Double_t GP_px2, GP_py2, GP_pz2; // Generated real photon from pi0 gen

    Double_t edep[1080];// energy deposition in each crystal
    Double_t X_sum, X_diff, X_BH; // cross sections

    t_g4->SetBranchAddress("evtNb", &evtNb); // Event number from event generator
    // Generated vertex position, the raster and beam offsets are included
    t_g4->SetBranchAddress("GV_x", &GV_x);
    t_g4->SetBranchAddress("GV_y", &GV_y);
    t_g4->SetBranchAddress("GV_z", &GV_z);
    // Generated Q2 and xB after the external energy loss correction before vertex
    t_g4->SetBranchAddress("GQ2", &GQ2);
    t_g4->SetBranchAddress("GxB", &GxB);
    // Generated t and phi from DVCS gen
    t_g4->SetBranchAddress("Gt", &Gt);
    t_g4->SetBranchAddress("Gphi", &Gphi);
    t_g4->SetBranchAddress("psf", &psf); // Phase Space Factor from DVCS gen
    // Beam energy of the kinematics
    t_g4->SetBranchAddress("RIE_px", &RIE_px);
    t_g4->SetBranchAddress("RIE_py", &RIE_py);
    t_g4->SetBranchAddress("RIE_pz", &RIE_pz);
    // Beam energy after the pre-vertex external radiation correction
    t_g4->SetBranchAddress("GIE_px", &GIE_px);
    t_g4->SetBranchAddress("GIE_py", &GIE_py);
    t_g4->SetBranchAddress("GIE_pz", &GIE_pz);
    // Generated scattered electron with the internal radiation correction
    t_g4->SetBranchAddress("GSE_px", &GSE_px);
    t_g4->SetBranchAddress("GSE_py", &GSE_py);
    t_g4->SetBranchAddress("GSE_pz", &GSE_pz);

    t_g4->SetBranchAddress("edep", edep); // Energy deposition in each crystal
    // Generated real photon from DVCS gen (not from Geant4)
    if(gen_type == 0){
        t_g4->SetBranchAddress("GP_px",  &GP_px);
        t_g4->SetBranchAddress("GP_py",  &GP_py);
        t_g4->SetBranchAddress("GP_pz",  &GP_pz);
    }
    else if(gen_type == 1){
        t_g4->SetBranchAddress("GP_px1",  &GP_px1);
        t_g4->SetBranchAddress("GP_py1",  &GP_py1);
        t_g4->SetBranchAddress("GP_pz1",  &GP_pz1);

        t_g4->SetBranchAddress("GP_px2",  &GP_px2);
        t_g4->SetBranchAddress("GP_py2",  &GP_py2);
        t_g4->SetBranchAddress("GP_pz2",  &GP_pz2);
    }

    t_g4->SetBranchAddress("X_sum", &X_sum); // XSecSum(0) from DVCS gen
    t_g4->SetBranchAddress("X_diff", &X_diff); // XSecDif() from DVCS gen, it will be 0 for pi0
    t_g4->SetBranchAddress("X_BH", &X_BH); // XSecSum(1) from DVCS gen, it will be 0 for pi0

    // From HMS simulation________________________________________________________________
    // Input file from the Geant4 simulation
    TString hms_inputDir = "worksim";
    TString hms_filename = Form("%s/hms_%s_%s_%d.root", hms_inputDir.Data(), Tar.Data(), kinc_param.c_str(), i_job);

    TFile *infile_hms = TFile::Open(Form("%s", hms_filename.Data()));
    if(!infile_hms){
        cout<<"Error: can not open "<<hms_filename.Data()<<endl;
        return;
    }
    TTree *t_hms = (TTree*)infile_hms->Get("h10");

    Double_t eventnumber; // Event number from HMS MC simulation
    // HMS variables from HMS MC simulation, only those with hms_stop_id = 0 (passed the HMS) are useful
    Double_t hms_stop_id; // 0: passed the full HMS; others: stopped in the simulation
    // Focal plane variables
    Double_t hsxfp;
    Double_t hsxpfp;
    Double_t hsyfp;
    Double_t hsypfp;
    // Variables at the target
    Double_t hsxptar;
    Double_t hsytar;
    Double_t hsyptar;
    Double_t hsdelta;

    t_hms->SetBranchAddress("eventnumber", &eventnumber); // Event number from HMS MC simulation
    t_hms->SetBranchAddress("hms_stop_id", &hms_stop_id); // Variables below have meaningful value only when hms_stop_id=0
    // Reconstructed scattered electon from HMS MC simulation
    t_hms->SetBranchAddress("MC_RSE_px", &RSE_px);
    t_hms->SetBranchAddress("MC_RSE_py", &RSE_py);
    t_hms->SetBranchAddress("MC_RSE_pz", &RSE_pz);
    // HMS coordinate variables
    t_hms->SetBranchAddress("hsxfp", &hsxfp);
    t_hms->SetBranchAddress("hsxpfp", &hsxpfp);
    t_hms->SetBranchAddress("hsyfp", &hsyfp);
    t_hms->SetBranchAddress("hsypfp", &hsypfp);
    t_hms->SetBranchAddress("hsxptar", &hsxptar);
    t_hms->SetBranchAddress("hsytar", &hsytar);
    t_hms->SetBranchAddress("hsyptar", &hsyptar);
    t_hms->SetBranchAddress("hsdelta", &hsdelta);


    // check if the number of events match
    if(t_g4->GetEntries() != t_hms->GetEntries()){
        cerr << "Error: The number of events in the Geant4 simulation and the HMS simulation do not match!" << endl;
        return;
    }

    // Output file for the reconstructed events________________________________________________________________
    TString outputDir = "rootfiles";
    TString output_filename = Form("%s/nps_hms_simulation_%s_%s_%d.root", outputDir.Data(), kinc_param.c_str(), Tar.Data(), i_job);

    TFile *output = new TFile(Form("%s", output_filename.Data()), "recreate");
    TTree *MC_dvcs = new TTree("MC_dvcs","TTree for DVCS simulation data"); // Create a new tree for reconstructed events

    // New branches for reconstructed photon and NPS clusters
    Double_t RV_z; // Reconstructed vertex z position
    Double_t clust_x, clust_y, clust_ene; // DVCS
    Int_t clust_size; // DVCS
    Double_t phot_px, phot_py, phot_pz; // DVCS

    Double_t clust_x1, clust_y1, clust_ene1; // pi0
    Double_t clust_x2, clust_y2, clust_ene2; // pi0
    Int_t clust_size1, clust_size2; // pi0
    Double_t phot_px1, phot_py1, phot_pz1; // pi0
    Double_t phot_px2, phot_py2, phot_pz2; // pi0

    Double_t Mx2, M;
    Double_t RQ2, Rphi, Rt, RxB;

    MC_dvcs->Branch("evtNb", &evtNb, "Event_Number/I");
    MC_dvcs->Branch("edep", edep, "Deposited_energy[1080]/D");
    MC_dvcs->Branch("GV_x", &GV_x, "Vertex_position_X_from_DVCS_gen/D");
    MC_dvcs->Branch("GV_y", &GV_y, "Vertex_position_Y_from_DVCS_gen/D");
    MC_dvcs->Branch("GV_z", &GV_z, "Vertex_position_Z_from_DVCS_gen/D");
    MC_dvcs->Branch("RIE_px", &RIE_px, "Initial_electron_Px_from_Geant4_is_beam_energy/D");
    MC_dvcs->Branch("RIE_py", &RIE_py, "Initial_electron_Py_from_Geant4_is_beam_energy/D");
    MC_dvcs->Branch("RIE_pz", &RIE_pz, "Initial_electron_Pz_from_Geant4_is_beam_energy/D");
    MC_dvcs->Branch("GIE_px", &GIE_px, "Initial_electron_Px_from_DVCS_gen/D");
    MC_dvcs->Branch("GIE_py", &GIE_py, "Initial_electron_Py_from_DVCS_gen/D");
    MC_dvcs->Branch("GIE_pz", &GIE_pz, "Initial_electron_Pz_from_DVCS_gen/D");
    MC_dvcs->Branch("GSE_px", &GSE_px, "Scattered_electron_Px_from_DVCS_gen/D");
    MC_dvcs->Branch("GSE_py", &GSE_py, "Scattered_electron_Py_from_DVCS_gen/D");
    MC_dvcs->Branch("GSE_pz", &GSE_pz, "Scattered_electron_Pz_from_DVCS_gen/D");
    if(gen_type == 0){
        MC_dvcs->Branch("GPh_px", &GP_px, "Real_photon_Px_from_DVCS_gen/D");
        MC_dvcs->Branch("GPh_py", &GP_py, "Real_photon_Py_from_DVCS_gen/D");
        MC_dvcs->Branch("GPh_pz", &GP_pz, "Real_photon_Pz_from_DVCS_gen/D");
    }
    else if(gen_type == 1){
        MC_dvcs->Branch("GPh_px1", &GP_px1, "Real_photon_Px_from_pi0_gen/D");
        MC_dvcs->Branch("GPh_py1", &GP_py1, "Real_photon_Py_from_pi0_gen/D");
        MC_dvcs->Branch("GPh_pz1", &GP_pz1, "Real_photon_Pz_from_pi0_gen/D");

        MC_dvcs->Branch("GPh_px2", &GP_px2, "Real_photon_Px_from_pi0_gen/D");
        MC_dvcs->Branch("GPh_py2", &GP_py2, "Real_photon_Py_from_pi0_gen/D");
        MC_dvcs->Branch("GPh_pz2", &GP_pz2, "Real_photon_Pz_from_pi0_gen/D");
    }
    
    MC_dvcs->Branch("GQ2", &GQ2, "Q2_from_DVCS_gen/D");
    MC_dvcs->Branch("GxB", &GxB, "xB_from_DVCS_gen/D");
    MC_dvcs->Branch("Gt", &Gt, "t_from_DVCS_gen/D");
    MC_dvcs->Branch("Gphi", &Gphi, "phi_from_DVCS_gen/D");
    MC_dvcs->Branch("psf", &psf, "Phase_Space_Factor_from_DVCS_gen/D");

    if(gen_type == 0){
        MC_dvcs->Branch("X_sum", &X_sum, "XSecSum(0)_from_DVCS_gen/D");
        MC_dvcs->Branch("X_diff", &X_diff, "XSecDif()_from_DVCS_gen/D");
        MC_dvcs->Branch("X_BH", &X_BH, "XSecSum(1)_from_DVCS_gen/D");
    }

    else if(gen_type == 1){
        MC_dvcs->Branch("X_sum", &X_sum, "XSecSum(0)_from_pi0_gen_(to_be_checked)/D");
    }

    MC_dvcs->Branch("hms_stop_id", &hms_stop_id, "HMS_stop_id/D");
    MC_dvcs->Branch("RV_z", &RV_z, "Reconstructed_vertex_Z_position_from_HMS_simulation/D");
    MC_dvcs->Branch("RSE_px", &RSE_px, "Scattered_electron_Px_from_in_HMS_simulation/D");
    MC_dvcs->Branch("RSE_py", &RSE_py, "Scattered_electron_Py_from_in_HMS_simulation/D");
    MC_dvcs->Branch("RSE_pz", &RSE_pz, "Scattered_electron_Pz_from_in_HMS_simulation/D");
    if(gen_type == 0){
        MC_dvcs->Branch("clust_x", &clust_x, "Cluster_X_position/D");
        MC_dvcs->Branch("clust_y", &clust_y, "Cluster_Y_position/D");
        MC_dvcs->Branch("clust_z", &caloDist, "Cluster_Z_position_is_the_NPS_distancen/D");
        MC_dvcs->Branch("clust_ene", &clust_ene, "Cluster_energy/D");
        MC_dvcs->Branch("clust_size", &clust_size, "Cluster_size/I");
    }
    else if(gen_type == 1){
        MC_dvcs->Branch("clust_x1", &clust_x1, "Cluster_X_position/D");
        MC_dvcs->Branch("clust_y1", &clust_y1, "Cluster_Y_position/D");
        MC_dvcs->Branch("clust_z1", &caloDist, "Cluster_Z_position_is_the_NPS_distancen/D");
        MC_dvcs->Branch("clust_ene1", &clust_ene1, "Cluster_energy/D");
        MC_dvcs->Branch("clust_size1", &clust_size1, "Cluster_size/I");

        MC_dvcs->Branch("clust_x2", &clust_x2, "Cluster_X_position/D");
        MC_dvcs->Branch("clust_y2", &clust_y2, "Cluster_Y_position/D");
        MC_dvcs->Branch("clust_z2", &caloDist, "Cluster_Z_position_is_the_NPS_distancen/D");
        MC_dvcs->Branch("clust_ene2", &clust_ene2, "Cluster_energy/D");
        MC_dvcs->Branch("clust_size2", &clust_size2, "Cluster_size/I");
    }
    MC_dvcs->Branch("M", &M, "Reconstructed_invariant_mass_should_be_0_if_only_one_photon/D");
    MC_dvcs->Branch("Mx2", &Mx2, "Missing_mass_squared_of_the_system/D");

    if(gen_type == 0){
        MC_dvcs->Branch("RPh_px", &phot_px, "Reconstructed_photon_Px/D");
        MC_dvcs->Branch("RPh_py", &phot_py, "Reconstructed_photon_Py/D");
        MC_dvcs->Branch("RPh_pz", &phot_pz, "Reconstructed_photon_Pz/D");
    }
    else if(gen_type == 1){
        MC_dvcs->Branch("RPh_px1", &phot_px1, "Reconstructed_photon1_Px/D");
        MC_dvcs->Branch("RPh_py1", &phot_py1, "Reconstructed_photon1_Py/D");
        MC_dvcs->Branch("RPh_pz1", &phot_pz1, "Reconstructed_photon1_Pz/D");

        MC_dvcs->Branch("RPh_px2", &phot_px2, "Reconstructed_photon2_Px/D");
        MC_dvcs->Branch("RPh_py2", &phot_py2, "Reconstructed_photon2_Py/D");
        MC_dvcs->Branch("RPh_pz2", &phot_pz2, "Reconstructed_photon2_Pz/D");
    }

    MC_dvcs->Branch("RQ2", &RQ2, "Reconstructed_Q2/D");
    MC_dvcs->Branch("RxB", &RxB, "Reconstructed_xB/D");
    MC_dvcs->Branch("Rt", &Rt, "Reconstructed_t/D");
    MC_dvcs->Branch("Rphi", &Rphi, "Reconstructed_phi/D");
    MC_dvcs->Branch("hsxfp", &hsxfp, "Focal_plane_X_from_HMS_simulation/D");
    MC_dvcs->Branch("hsxpfp", &hsxpfp, "Focal_plane_X'_from_HMS_simulation/D");
    MC_dvcs->Branch("hsyfp", &hsyfp, "Focal_plane_Y_from_HMS_simulation/D");
    MC_dvcs->Branch("hsypfp", &hsypfp, "Focal_plane_Y'_from_HMS_simulation/D");
    MC_dvcs->Branch("hsxptar", &hsxptar, "X_at_target_from_HMS_simulation/D");
    MC_dvcs->Branch("hsytar", &hsytar, "Y_at_target_from_HMS_simulation/D");
    MC_dvcs->Branch("hsyptar", &hsyptar, "Y'_at_target_from_HMS_simulation/D");
    MC_dvcs->Branch("hsdelta", &hsdelta, "Delta_at_target_from_HMS_simulation/D");

    //w0 and shower depth information just in case
    Float_t a, x, y, x_corr, y_corr; // For DVCS
    Float_t a1, x1, y1, x_corr1, y_corr1; // For pi0 photon1
    Float_t a2, x2, y2, x_corr2, y_corr2; // For pi0 photon2
    if(gen_type == 0){
        MC_dvcs->Branch("a", &a, "shower_depth_along_photon_trajectory/F");
        MC_dvcs->Branch("x_corr", &x_corr, "cluster x after shower depth correction/F");
        MC_dvcs->Branch("y_corr", &y_corr, "cluster y after shower depth correction/F");
    }
    else if(gen_type == 1){
        MC_dvcs->Branch("a1", &a1, "shower_depth_along_photon_trajectory/F");
        MC_dvcs->Branch("x_corr1", &x_corr1, "cluster x after shower depth correction/F");
        MC_dvcs->Branch("y_corr1", &y_corr1, "cluster y after shower depth correction/F");

        MC_dvcs->Branch("a2", &a2, "shower_depth_along_photon_trajectory/F");
        MC_dvcs->Branch("x_corr2", &x_corr2, "cluster x after shower depth correction/F");
        MC_dvcs->Branch("y_corr2", &y_corr2, "cluster y after shower depth correction/F");
    }
    
    // Start reconstruction and analysis loop____________________________________________________________
    TDVCSEvent *dvcs_evt = new TDVCSEvent(run_number[0]);
    TCaloEvent *calo_evt = new TCaloEvent(run_number[0]);
    TLorentzVector *L_calo_phot;// for DVCS
    TLorentzVector *L_calo_phot1, *L_calo_phot2;// for pi0

    Double_t CountEvt = 0; // for the reconstuction with different dead block configurations
    Int_t NPSconfig = 0; // NPS dead block configuration index

    Int_t nevt = t_hms->GetEntries();
    // Int_t nevt = 1000000; // for test
    for(int ievt = 0; ievt < nevt; ievt++){
        if(ievt%1000 == 0) cout<<ievt<<"/"<<nevt<<endl;
        t_g4->GetEntry(ievt);
        t_hms->GetEntry(ievt);
        if(evtNb != eventnumber) {
            cerr << "Error: Event number mismatch!!!!!!!!!!" << ievt<< endl;
            return;
        }
        // initialize the reconstructed variables, updated only when hms_stop_id = 0 and Nb_clust >= 1
        RV_z = -999;
        
        clust_x = -999; clust_y = -999; clust_ene = -999; clust_size = -999;// DVCS
        phot_px = -999; phot_py = -999; phot_pz = -999;// DVCS
         
        clust_x1 = -999; clust_y1 = -999; clust_ene1 = -999; clust_size1 = -999;// pi0
        phot_px1 = -999; phot_py1 = -999; phot_pz1 = -999;// pi0
        clust_x2 = -999; clust_y2 = -999; clust_ene2 = -999; clust_size2 = -999;// pi0
        phot_px2 = -999; phot_py2 = -999; phot_pz2 = -999;// pi0

        M = -999; Mx2 = -999;// DVCS and pi0
        RQ2 = -999; Rphi = -999; Rt = -999; RxB = -999;// DVCS and pi0

        a = -999; x_corr = -999; y_corr = -999;// DVCS
        a1 = -999; x_corr1 = -999; y_corr1 = -999;// pi0
        a2 = -999; x_corr2 = -999; y_corr2 = -999;// pi0

        // For different dead block configurations, assign the caloMaskBlock array accordingly
        CountEvt++;
        // cout<<CountEvt/nevt<<", "<<total_charge_frac[NPSconfig]<<endl;
        if(CountEvt/nevt > total_charge_frac[NPSconfig] && NPSconfig < run_number.size()-1){
            NPSconfig++;
            caloMaskBlock_temp = db->GetEntry_i("CALO_flag_MaskBlock", run_number[NPSconfig]);
            // Convert the mask block information to the simulation numbering scheme
            for(Int_t i = 0; i < 1080; i++){
                caloMaskBlock[i] = caloMaskBlock_temp[bnConv_OldToNew(i)];
            }
            cout<<"Now at event "<<ievt<<", switch to run "<<run_number[NPSconfig]<<" for dead block configuration #"<<NPSconfig+1<<endl;
            cout<<"========== NPS mask blocks =========="<<endl;
            for(int i = 0; i < 1080; i++) if(caloMaskBlock_temp[i] == 1) cout<<i<<" "; // output in NPS numbering scheme
            cout<<"======================================"<<endl;
            cout<<endl;
        }

        if(hms_stop_id == 0){ // Reconstruct only when the event passed the HMS simulation
            RV_z = GetVertexZ(GV_x, HMS_angle, hsytar, hsyptar); // reconstruction of vertex z

            dvcs_evt->SetVertex(GV_x, GV_y, RV_z);//[cm]
            dvcs_evt->GetGeometry()->SetCaloTheta(-1*caloAngle);
            dvcs_evt->GetGeometry()->SetCaloDist(caloDist);
            
            L_calo_phot = new TLorentzVector(); // DVCS
            L_calo_phot1 = new TLorentzVector(); // pi0
            L_calo_phot2 = new TLorentzVector(); // pi0

            dvcs_evt->SetCaloEvent(calo_evt);

            // cout<<"Event "<<ievt<<": edep in each crystal (MeV): ";
            // for(int i = 0; i < 1080; i++) cout<<edep[i]*(1e-3)<<" ";
            // cout<<endl;
            // cout<<"================================================"<<endl;
            // continue;
            
            // nRecEvt++; // used for different dead block configurations
            for (int i = 0 ; i < 1080 ; i++){ // NOTE!! This should be in simulation numbering scheme!!!!!!!!!!!!!    
                TCaloBlock* block=calo_evt->AddBlock(i);

                //energy (GeV) and time (ns) (currently, no time info in Geant4 simulation. default as 0 for now.)
                if(caloMaskBlock[i]) block->AddPulse(0, 0.); // No pulse in dead blocks
                else block->AddPulse(1.04*edep[i]*(1e-3), 0.); // Energy after calibration for simulation
            }

            calo_evt->TriggerSim(clusThr); // set higher if too many blocks in the event
            calo_evt->DoClustering(-3., 3.);//-3ns to +3ns windows
            Int_t Nb_clust = calo_evt->GetNbClusters();
            // cout<<"Number of clusters: "<<Nb_clust<<endl;
            for(int i = 0 ; i < Nb_clust ; i++ ){
                TCaloCluster *clus = calo_evt->GetCluster(i);
                clus->Analyze(1., -3., 3.);//"true||false", "time_min(ns)", "time_max(ns)", "weight"
            }

            if(gen_type == 0 && Nb_clust == 1) { // Only process events that have 1 clusters in NPS for DVCS
                clust_x   = calo_evt->GetCluster(0)->GetX();
                clust_y   = calo_evt->GetCluster(0)->GetY();
                clust_ene = calo_evt->GetCluster(0)->GetE();
                clust_size   = calo_evt->GetCluster(0)->GetClusSize();
                
                // cout<<calo_evt->GetCluster(0)->GetClusSize()<<endl;
                // cout<<"Cluster (x,y): ("<<clust_x<<", "<<clust_y<<") with energy "<<clust_ene<<" GeV and size "<<clust_size<<endl;
                *L_calo_phot = dvcs_evt->GetPhoton(0, 7, 0, a, x_corr, y_corr);
                phot_px = L_calo_phot->Px();//GeV
                phot_py = L_calo_phot->Py();
                phot_pz = L_calo_phot->Pz();
                // cout<<"photon: ("<<phot_px<<", "<<phot_py<<", "<<phot_pz<<")"<<endl;
                // if(clust_ene>0) cout<<a<<endl;
                TVector3 k(0, 0, RIE_pz);
                TVector3 kp(RSE_px, RSE_py, RSE_pz);
                TVector3 qp(phot_px, phot_py, phot_pz);

                RQ2 = 2*k.Mag()*kp.Mag()*(1 - cos(k.Angle(kp)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))

                TVector3 q = k - kp;

                TVector3 v1=q.Cross(kp);
                TVector3 v2=q.Cross(qp);
                Rphi=v1.Angle(v2);
                if(q.Dot(v1.Cross(v2))<0) Rphi=2.*TMath::Pi()-Rphi;

                //debug
                // cout<<"k: ("<<k.Px()<<", "<<k.Py()<<", "<<k.Pz()<<")"<<endl;
                // cout<<"kp: ("<<kp.Px()<<", "<<kp.Py()<<", "<<kp.Pz()<<")"<<endl;
                // cout<<"qp: ("<<qp.Px()<<", "<<qp.Py()<<", "<<qp.Pz()<<")"<<endl;
                // cout<<"q: ("<<q.Px()<<", "<<q.Py()<<", "<<q.Pz()<<")"<<endl;
                // cout<<"q cross kp: ("<<v1.Px()<<", "<<v1.Py()<<", "<<v1.Pz()<<")"<<endl;
                // cout<<"q cross qp: ("<<v2.Px()<<", "<<v2.Py()<<", "<<v2.Pz()<<")"<<endl;

                Double_t cos=(q.Dot(qp))/(q.Mag()*qp.Mag());
                Double_t nu=RIE_pz-kp.Mag();
                Rt=(RQ2*M_targ+2.*nu*M_targ*(nu-sqrt(nu*nu+RQ2)*cos))/(sqrt(nu*nu+RQ2)*cos-nu-M_targ);
                RxB = RQ2/(2*M_targ*nu);

                TLorentzVector Lb, Lp, Lgp, Lkp, L_mis;
                Double_t RSE_E = sqrt(RSE_px*RSE_px+RSE_py*RSE_py+RSE_pz*RSE_pz+0.000511*0.000511); // scattered electron energy from HMS simulation
                Lb.SetPxPyPzE(0, 0, RIE_pz, RIE_pz);
                Lp.SetPxPyPzE(0, 0, 0, M_targ);
                Lkp.SetPxPyPzE(RSE_px, RSE_py, RSE_pz, RSE_E);
                Lgp = *L_calo_phot;

                Mx2 = (Lb + Lp - Lgp - Lkp).M2();
                M = Lgp.M();

                // cout<<"Lb: ("<<Lb.Px()<<", "<<Lb.Py()<<", "<<Lb.Pz()<<", "<<Lb.E()<<")"<<endl;
                // cout<<"Lp: ("<<Lp.Px()<<", "<<Lp.Py()<<", "<<Lp.Pz()<<", "<<Lp.E()<<")"<<endl;
                // cout<<"Lkp: ("<<Lkp.Px()<<", "<<Lkp.Py()<<", "<<Lkp.Pz()<<", "<<Lkp.E()<<")"<<endl;
                // cout<<"Lgp: ("<<Lgp.Px()<<", "<<Lgp.Py()<<", "<<Lgp.Pz()<<", "<<Lgp.E()<<")"<<endl;
                // cout<<"Mx2: "<<Mx2<<endl;
            }//Nb_clus>=1 && DVCS

            if(gen_type == 1) { // pi0
                if(Nb_clust == 1){ // 1 cluster event
                    clust_x1   = calo_evt->GetCluster(0)->GetX();
                    clust_y1   = calo_evt->GetCluster(0)->GetY();
                    clust_ene1 = calo_evt->GetCluster(0)->GetE();
                    clust_size1   = calo_evt->GetCluster(0)->GetClusSize();
                    
                    // cout<<calo_evt->GetCluster(0)->GetClusSize()<<endl;
                    // cout<<"Cluster (x,y): ("<<clust_x1<<", "<<clust_y1<<") with energy "<<clust_ene1<<" GeV and size "<<clust_size1<<endl;
                    *L_calo_phot1 = dvcs_evt->GetPhoton(0, 7, 0, a1, x_corr1, y_corr1);
                    phot_px1 = L_calo_phot1->Px();//GeV
                    phot_py1 = L_calo_phot1->Py();
                    phot_pz1 = L_calo_phot1->Pz();
                    // cout<<"photon: ("<<phot_px1<<", "<<phot_py1<<", "<<phot_pz1<<")"<<endl;
                    // if(clust_ene1>0) cout<<a<<endl;
                    TVector3 k(0, 0, RIE_pz);
                    TVector3 kp(RSE_px, RSE_py, RSE_pz);
                    TVector3 qp(phot_px1, phot_py1, phot_pz1);

                    RQ2 = 2*k.Mag()*kp.Mag()*(1 - cos(k.Angle(kp)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))

                    TVector3 q = k - kp;

                    TVector3 v1=q.Cross(kp);
                    TVector3 v2=q.Cross(qp);
                    Rphi=v1.Angle(v2);
                    if(q.Dot(v1.Cross(v2))<0) Rphi=2.*TMath::Pi()-Rphi;

                    //debug
                    // cout<<"k: ("<<k.Px()<<", "<<k.Py()<<", "<<k.Pz()<<")"<<endl;
                    // cout<<"kp: ("<<kp.Px()<<", "<<kp.Py()<<", "<<kp.Pz()<<")"<<endl;
                    // cout<<"qp: ("<<qp.Px()<<", "<<qp.Py()<<", "<<qp.Pz()<<")"<<endl;
                    // cout<<"q: ("<<q.Px()<<", "<<q.Py()<<", "<<q.Pz()<<")"<<endl;
                    // cout<<"q cross kp: ("<<v1.Px()<<", "<<v1.Py()<<", "<<v1.Pz()<<")"<<endl;
                    // cout<<"q cross qp: ("<<v2.Px()<<", "<<v2.Py()<<", "<<v2.Pz()<<")"<<endl;

                    Double_t cos=(q.Dot(qp))/(q.Mag()*qp.Mag());
                    Double_t nu=RIE_pz-kp.Mag();
                    Rt=(RQ2*M_targ+2.*nu*M_targ*(nu-sqrt(nu*nu+RQ2)*cos))/(sqrt(nu*nu+RQ2)*cos-nu-M_targ);
                    RxB = RQ2/(2*M_targ*nu);

                    TLorentzVector Lb, Lp, Lgp, Lkp, L_mis;
                    Double_t RSE_E = sqrt(RSE_px*RSE_px+RSE_py*RSE_py+RSE_pz*RSE_pz+0.000511*0.000511); // scattered electron energy from HMS simulation
                    Lb.SetPxPyPzE(0, 0, RIE_pz, RIE_pz);
                    Lp.SetPxPyPzE(0, 0, 0, M_targ);
                    Lkp.SetPxPyPzE(RSE_px, RSE_py, RSE_pz, RSE_E);
                    Lgp = *L_calo_phot1;

                    Mx2 = (Lb + Lp - Lgp - Lkp).M2();
                    M = Lgp.M();

                    // cout<<"Lb: ("<<Lb.Px()<<", "<<Lb.Py()<<", "<<Lb.Pz()<<", "<<Lb.E()<<")"<<endl;
                    // cout<<"Lp: ("<<Lp.Px()<<", "<<Lp.Py()<<", "<<Lp.Pz()<<", "<<Lp.E()<<")"<<endl;
                    // cout<<"Lkp: ("<<Lkp.Px()<<", "<<Lkp.Py()<<", "<<Lkp.Pz()<<", "<<Lkp.E()<<")"<<endl;
                    // cout<<"Lgp: ("<<Lgp.Px()<<", "<<Lgp.Py()<<", "<<Lgp.Pz()<<", "<<Lgp.E()<<")"<<endl;
                    // cout<<"Mx2: "<<Mx2<<endl;
                } // 1 cluster

                if(Nb_clust == 2){ // 2 cluster event
                    clust_x1   = calo_evt->GetCluster(0)->GetX();
                    clust_y1   = calo_evt->GetCluster(0)->GetY();
                    clust_ene1 = calo_evt->GetCluster(0)->GetE();
                    clust_size1   = calo_evt->GetCluster(0)->GetClusSize();

                    clust_x2   = calo_evt->GetCluster(1)->GetX();
                    clust_y2   = calo_evt->GetCluster(1)->GetY();
                    clust_ene2 = calo_evt->GetCluster(1)->GetE();
                    clust_size2   = calo_evt->GetCluster(1)->GetClusSize();
                    
                    // cout<<calo_evt->GetCluster(0)->GetClusSize()<<endl;
                    // cout<<"Cluster1 (x,y): ("<<clust_x1<<", "<<clust_y1<<") with energy "<<clust_ene1<<" GeV and size "<<clust_size1<<endl;
                    // cout<<"Cluster2 (x,y): ("<<clust_x2<<", "<<clust_y2<<") with energy "<<clust_ene2<<" GeV and size "<<clust_size2<<endl;
                    *L_calo_phot1 = dvcs_evt->GetPhoton(0, 7, 0, a1, x_corr1, y_corr1);
                    phot_px1 = L_calo_phot1->Px();//GeV
                    phot_py1 = L_calo_phot1->Py();
                    phot_pz1 = L_calo_phot1->Pz();

                    *L_calo_phot2 = dvcs_evt->GetPhoton(1, 7, 0, a2, x_corr2, y_corr2);
                    phot_px2 = L_calo_phot2->Px();//GeV
                    phot_py2 = L_calo_phot2->Py();
                    phot_pz2 = L_calo_phot2->Pz();
                    // cout<<"photon: ("<<phot_px1<<", "<<phot_py1<<", "<<phot_pz1<<")"<<endl;
                    // cout<<"photon: ("<<phot_px2<<", "<<phot_py2<<", "<<phot_pz2<<")"<<endl;
                    // if(clust_ene1>0) cout<<a<<endl;
                    TVector3 k(0, 0, RIE_pz);
                    TVector3 kp(RSE_px, RSE_py, RSE_pz);
                    TVector3 qp(phot_px1+phot_px2, phot_py1+phot_py2, phot_pz1+phot_pz2); // pi0 momenta for t and phi

                    RQ2 = 2*k.Mag()*kp.Mag()*(1 - cos(k.Angle(kp)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))

                    TVector3 q = k - kp;

                    TVector3 v1=q.Cross(kp);
                    TVector3 v2=q.Cross(qp);
                    Rphi=v1.Angle(v2);
                    if(q.Dot(v1.Cross(v2))<0) Rphi=2.*TMath::Pi()-Rphi;

                    //debug
                    // cout<<"k: ("<<k.Px()<<", "<<k.Py()<<", "<<k.Pz()<<")"<<endl;
                    // cout<<"kp: ("<<kp.Px()<<", "<<kp.Py()<<", "<<kp.Pz()<<")"<<endl;
                    // cout<<"qp: ("<<qp.Px()<<", "<<qp.Py()<<", "<<qp.Pz()<<")"<<endl;
                    // cout<<"q: ("<<q.Px()<<", "<<q.Py()<<", "<<q.Pz()<<")"<<endl;
                    // cout<<"q cross kp: ("<<v1.Px()<<", "<<v1.Py()<<", "<<v1.Pz()<<")"<<endl;
                    // cout<<"q cross qp: ("<<v2.Px()<<", "<<v2.Py()<<", "<<v2.Pz()<<")"<<endl;

                    Double_t cos=(q.Dot(qp))/(q.Mag()*qp.Mag());
                    Double_t nu=RIE_pz-kp.Mag();
                    Rt=(RQ2*M_targ+2.*nu*M_targ*(nu-sqrt(nu*nu+RQ2)*cos))/(sqrt(nu*nu+RQ2)*cos-nu-M_targ);
                    RxB = RQ2/(2*M_targ*nu);

                    TLorentzVector Lb, Lp, Lgp, Lkp, L_mis;
                    Double_t RSE_E = sqrt(RSE_px*RSE_px+RSE_py*RSE_py+RSE_pz*RSE_pz+0.000511*0.000511); // scattered electron energy from HMS simulation
                    Lb.SetPxPyPzE(0, 0, RIE_pz, RIE_pz);
                    Lp.SetPxPyPzE(0, 0, 0, M_targ);
                    Lkp.SetPxPyPzE(RSE_px, RSE_py, RSE_pz, RSE_E);
                    Lgp = *L_calo_phot1+*L_calo_phot2;

                    Mx2 = (Lb + Lp - Lgp - Lkp).M2();
                    M = Lgp.M();
                    // cout<<"Lb: ("<<Lb.Px()<<", "<<Lb.Py()<<", "<<Lb.Pz()<<", "<<Lb.E()<<")"<<endl;
                    // cout<<"Lp: ("<<Lp.Px()<<", "<<Lp.Py()<<", "<<Lp.Pz()<<", "<<Lp.E()<<")"<<endl;
                    // cout<<"Lkp: ("<<Lkp.Px()<<", "<<Lkp.Py()<<", "<<Lkp.Pz()<<", "<<Lkp.E()<<")"<<endl;
                    // cout<<"Lgp: ("<<Lgp.Px()<<", "<<Lgp.Py()<<", "<<Lgp.Pz()<<", "<<Lgp.E()<<")"<<endl;
                    // cout<<"Mx2: "<<Mx2<<endl;
                }
            }//pi0

            calo_evt->Reset();
        } //hms_stop_id=0

        MC_dvcs->Fill();
    }
    output->cd();
    MC_dvcs->Write();
    output->Close();
}

void GetRunRange(const string &kinc_param, vector<int> &minrun, vector<int> &maxrun) {
    minrun.clear();
    maxrun.clear();

    // Get the kinematic list from mc-single-arm package
    const char* folder_hms_c = gSystem->Getenv("folder_hms");
    if(!folder_hms_c){
        cerr << "Error: environment variable folder_hms not set!" << endl;
        return;
    }
    const string path_folder_hms(folder_hms_c);
    const string filename = path_folder_hms + "/KinC_list.txt";
    ifstream infile(filename);

    if (!infile.is_open()) {
        cerr << "Error: can't not open " << filename << endl;
        return;
    }

    // Get the list of begin and end run, put into the vectors
    string line;
    if (!getline(infile, line)) return;

    string last_kinc_name = "";
    while (getline(infile, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        string token;
        
        char delimiter = '\t'; 

        getline(ss, token, delimiter); 
        string current_kinc = token;
        
        if (current_kinc.find_first_not_of(" \t") != string::npos) {
            current_kinc = current_kinc.substr(current_kinc.find_first_not_of(" \t"));
        } else {
            current_kinc = "";
        }
        
        string begin_run_str;
        getline(ss, begin_run_str, delimiter);

        string end_run_str;
        getline(ss, end_run_str, delimiter);
        
        string target_kinc_name;
        if (!current_kinc.empty()) {
            last_kinc_name = current_kinc;
            target_kinc_name = current_kinc;
        } else {
            target_kinc_name = last_kinc_name;
        }

        if (target_kinc_name == kinc_param) {
            try {
                int begin_run = stoi(begin_run_str);
                int end_run = stoi(end_run_str);

                minrun.push_back(begin_run);
                maxrun.push_back(end_run);
            } catch (const exception& e) {
                cerr << "Warning: failed transfer from string to int, skip this line: " << line << std::endl;
            }
        }
    }

    infile.close();
}

// Get the run numbers and their total charge from good production runs
void GetRunNumbersAndCharge(int target_flag, vector<int> minrun, vector<int> maxrun, vector <int> &run_number, vector <double> &total_charge){
    Int_t Prod_Flag = 0;
    Int_t Quality_Flag = -2;
    Int_t Target_amu = 0;
    Double_t Charge = 0;

    Int_t nSubKine = minrun.size();
    if(nSubKine != maxrun.size()){
        cout << "Error: minrun and maxrun vectors have different sizes!" << endl;
        return;
    }

    for(int i = 0; i < nSubKine; i++){
        for(int run = minrun[i]; run <= maxrun[i]; run++){
            Prod_Flag = *db->GetEntry_i("RUN_flag_IsProdRun", run);
            Quality_Flag = *db->GetEntry_i("RUN_flag_Quality", run);
            Target_amu = *db->GetEntry_d("TARGET_param_Amu", run);
            Charge = *db->GetEntry_d("BEAM_param_TotalCharge", run);
            
            if(Prod_Flag == 1 && Quality_Flag == 0){
                if(target_flag == 0){// LH2
                    if(0 < Target_amu && Target_amu < 1.5){
                        run_number.push_back(run);
                        total_charge.push_back(Charge);
                    }
                }
                else if(target_flag == 1){// LD2
                    if(1.5 < Target_amu && Target_amu < 2.5){
                        run_number.push_back(run);
                        total_charge.push_back(Charge);
                    }
                }
            }
        }
    }
}

// Reconstruction of vertex z
double GetVertexZ(double GV_x, double HMS_Angle, double hsytar, double hsyptar){
    Double_t vertex_z;
    Double_t HMS_Angle_cos = cos(HMS_Angle);
    Double_t HMS_Angle_sin = sin(HMS_Angle);
    if(hsytar<-1e5) vertex_z = -1e6;
    else{
        vertex_z = hsytar + 0.1*(0.52-0.012*HMS_Angle+0.002*HMS_Angle*HMS_Angle) - GV_x*( HMS_Angle_cos + hsyptar*HMS_Angle_sin );
        // vertex_z = hsytar - GV_x*( HMS_Angle_cos + hsyptar*HMS_Angle_sin );
        vertex_z = vertex_z / ( HMS_Angle_sin - hsyptar*HMS_Angle_cos );
    }
    return vertex_z;
}

// Functions to convert the numbering scheme of NPS
Int_t bnConv_OldToNew(int ibn_old) // Convert the PMT number to the current version
{
  const Int_t ncol = 30; // number of columns
  const Int_t nrow = 36; // number of rows
  Int_t irow = ibn_old % nrow;
  Int_t icol = 29 - ((ibn_old - irow) / nrow);
  Int_t ibn_new = ncol * irow + icol;

  return ibn_new;
}

Int_t bnConv_NewToOld(int ibn_new) // Convert the PMT number to the simulation version
{
  const Int_t ncol = 30; // number of columns
  const Int_t nrow = 36; // number of rows
  Int_t icol = 29 - (ibn_new % ncol);
  Int_t irow = (ibn_new - ibn_new % ncol) / ncol;
  Int_t ibn_old = nrow * icol + irow;

  return ibn_old;
}