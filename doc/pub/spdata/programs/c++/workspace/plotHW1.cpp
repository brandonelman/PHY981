// Liquid drop model code
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "TLegend.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TCanvas.h"

//Calculate liquid drop model binding energies in MeV
//Uses formula:
//BE(A,Z) = C1*A - C2*A^(2,3) - C3*(Z^2/A^(1/3)) - C4(A-2Z)^2/A 
double calculateBindingEnergyLD(int Z, int A){
  const double VOL_CONST  = 15.49;
  const double SURF_CONST = 17.23;
  const double COUL_CONST = 0.697;
  const double SYMM_CONST =  22.6;

  double binding_energy = VOL_CONST * A;
  binding_energy -= SURF_CONST * pow(A, 2./3.);
  binding_energy -= COUL_CONST * pow(Z, 2.) / pow(A,1./3.); 
  binding_energy -= SYMM_CONST * pow((A-2*Z), 2.) / A;
  return binding_energy;
}

double calculateBindingEnergyLD(int Z, int A, const double VOL_CONST, const double SURF_CONST,
                                const double COUL_CONST, const double SYMM_CONST){
  double binding_energy = VOL_CONST * A;
  binding_energy -= SURF_CONST * pow(A, 2./3.);
  binding_energy -= COUL_CONST * pow(Z, 2.) / pow(A,1./3.); 
  binding_energy -= SYMM_CONST * pow((A-2*Z), 2.) / A;
  return binding_energy;
}

//This function will take as input the file bedata.dat and output the plots 
//necessary for Hw1.
//
//The input file is structured as: 
//Z  A  binding_energy error J N  (where function of J is unknown)
void plotHW1(std::string input_file_name, std::string output_file_name){
  
  int lines_read = 0;

  const int LINES_EXPECTED = 2375; //expected lines from wc -l
  const int Z_OXYGEN = 8;
  const int Z_CALCIUM = 20;
  const int Z_TIN = 50;
  const int Z_LEAD = 82;

  const int MAX_Z = 110;
  const int MAX_N = 159;
  const int MAX_A = 269;

  const double VOL_CONST = 15.49;
  const double SURF_CONST = 17.23;
  const double COUL_CONST = 0.697;
  const double SYMM_CONST = 22.6;

  int Z[LINES_EXPECTED]; //# of protons
  int A[LINES_EXPECTED]; //# of nucleons
  double binding_energy[LINES_EXPECTED];//in MeV
  double be_error[LINES_EXPECTED];
  int garb[LINES_EXPECTED];//unclear what this is; always equals  1
  int N[LINES_EXPECTED];//# of neutrons

  //Referenced from 0; so neutron_sep_energy[0][0] would correspond to the deuteron
  double neutron_sep_energy[MAX_Z][MAX_N];
  double neutron_sep_energy_LD[MAX_Z][MAX_N]; // Liquid drop version
  for (int i = 0; i < MAX_Z; i++){
    for (int j = 0; j < MAX_N; j++){
      neutron_sep_energy[i][j] = 0.0;
      neutron_sep_energy_LD[i][j] = 0.0;
    }
  }

  double binding_energy_per_nucleon[MAX_A];
  double binding_energy_per_nucleon_LD[MAX_A];
  double binding_energy_per_nucleon_LD_vol[MAX_A];
  double binding_energy_per_nucleon_LD_volsurf[MAX_A];
  double binding_energy_per_nucleon_LD_volsurfcoul[MAX_A];
  for (int i = 0; i < MAX_A; i++){
    binding_energy_per_nucleon[i] = 0.0;
    binding_energy_per_nucleon_LD[i] = 0.0;
    binding_energy_per_nucleon_LD_vol[i] = 0.0;
    binding_energy_per_nucleon_LD_volsurf[i] = 0.0;
    binding_energy_per_nucleon_LD_volsurfcoul[i] = 0.0;
  }
  
  int n_drip_line[120];
  int n_points[MAX_Z];
  TGraph *LD_sep_energy_plot[MAX_Z];
  TGraph *sep_energy_plot[MAX_Z];
  TMultiGraph *sep_energy_graphs = new TMultiGraph();
  TMultiGraph *be_graphs_LD = new TMultiGraph();
  TGraph *n_equals_z_line = new TGraph();
  TGraph *n_dripline_plot = new TGraph();

  TMultiGraph *dripline_plots = new TMultiGraph();

  for (int i = 0; i < MAX_Z; i++){
    sep_energy_plot[i] = new TGraph();
    LD_sep_energy_plot[i] = new TGraph();
    n_points[i] = 0;
    n_drip_line[i] = 0; 
    if (i < 120){
      n_equals_z_line->SetPoint(i, i+1, i+1);
    }
  }

  //For parsing input file
  std::string line;
  std::ifstream input_file;


  input_file.open(input_file_name.c_str());
  if (!input_file){
    std::cout << "Failed to open input file with name: " << input_file_name << std::endl;
    return;
  }
  
  //First lets parse the input file 
  while(std::getline(input_file,line)){
    std::stringstream ss(line);
    ss >> Z[lines_read];
    ss >> A[lines_read];
    ss >> binding_energy[lines_read];
    ss >> be_error[lines_read];
    ss >> garb[lines_read];
    ss >> N[lines_read];
    if (N[lines_read] + Z[lines_read] != A[lines_read]){
      std::cout << "This makes no sense! Proton Number + Neutron Number != mass number!" << std::endl;
    }
    lines_read++;
  }
  
  if (lines_read != LINES_EXPECTED){
    std::cout << "Read wrong numbers of lines: " << lines_read << " compared to expected: "
              << LINES_EXPECTED << std::endl;
    input_file.close();
    return;
  }
  input_file.close();
  
  int cur_line = 0;
  while (cur_line < LINES_EXPECTED){
    if (binding_energy_per_nucleon[A[cur_line]-1] < binding_energy[cur_line]/A[cur_line]){
      binding_energy_per_nucleon[A[cur_line]-1] = binding_energy[cur_line]/A[cur_line];
    }

    double binding_energy_LD_temp = calculateBindingEnergyLD(Z[cur_line],A[cur_line]) / A[cur_line];
    if (binding_energy_per_nucleon_LD[A[cur_line]-1] < binding_energy_LD_temp){
      binding_energy_per_nucleon_LD[A[cur_line]-1] = binding_energy_LD_temp;
    }

    double binding_energy_LD_temp_vol      = calculateBindingEnergyLD(Z[cur_line],A[cur_line], VOL_CONST, 0,0,0) / A[cur_line];
    if (binding_energy_per_nucleon_LD_vol[A[cur_line]-1] < binding_energy_LD_temp_vol){
      binding_energy_per_nucleon_LD_vol[A[cur_line]-1] = binding_energy_LD_temp_vol;
    }
    double binding_energy_LD_temp_volsurf     = calculateBindingEnergyLD(Z[cur_line],A[cur_line], VOL_CONST, SURF_CONST,0,0) / A[cur_line];
    if (binding_energy_per_nucleon_LD_volsurf[A[cur_line]-1] < binding_energy_LD_temp_volsurf){
      binding_energy_per_nucleon_LD_volsurf[A[cur_line]-1] = binding_energy_LD_temp_volsurf;
    }
    double binding_energy_LD_temp_volsurfcoul = calculateBindingEnergyLD(Z[cur_line],A[cur_line], VOL_CONST, SURF_CONST,COUL_CONST,0) / A[cur_line];
    if (binding_energy_per_nucleon_LD_volsurfcoul[A[cur_line]-1] < binding_energy_LD_temp_volsurfcoul){
      binding_energy_per_nucleon_LD_volsurfcoul[A[cur_line]-1] = binding_energy_LD_temp_volsurfcoul;
    }

    cur_line++;
    if (Z[cur_line] == Z[cur_line-1]){
      if(N[cur_line]-N[cur_line-1] == 1){
        neutron_sep_energy[Z[cur_line]-1][N[cur_line]-1] = binding_energy[cur_line]-binding_energy[cur_line-1];
        //For N 
        double binding_energy_LD1 = calculateBindingEnergyLD(Z[cur_line], A[cur_line]); 
        //For N-1
        double binding_energy_LD2 = calculateBindingEnergyLD(Z[cur_line-1], A[cur_line-1]);
        double separation_energy_LD = binding_energy_LD1 - binding_energy_LD2;
//      if (separation_energy_LD < 0 && n_drip_line[Z[cur_line]-1] == 0){
//        std::cout << "Success! n_drip_line for Z = "<< Z[cur_line] << " is N = " << N[cur_line] << std::endl;
//        n_drip_line[Z[cur_line]-1] = N[cur_line];
//      }
        neutron_sep_energy_LD[Z[cur_line]-1][N[cur_line]-1] = separation_energy_LD;
      }
    }
    continue;
  }
  

  bool found;
  for (int i = 1; i <= 120; i++){
    found = false;
    int num_protons = i;
    int num_neutrons = i;
    while (!found){
      double binding_energy_LD1  =  calculateBindingEnergyLD(num_protons, (num_protons+num_neutrons)); 
      double binding_energy_LD2  =  calculateBindingEnergyLD(num_protons, (num_protons+num_neutrons-1)); 
      double separation_energy_LD = binding_energy_LD1 - binding_energy_LD2;
      if (separation_energy_LD < 0){
        found = true;
        n_drip_line[num_protons] = num_neutrons;
        std::cout << "Drip line for Z = " << num_protons << " is N = " << n_drip_line[num_protons] << std::endl;
      }
      else{
        num_neutrons++;
      }
    }
    n_dripline_plot->SetPoint(i-1, n_drip_line[num_protons], num_protons);
  }
  
  for (int i = 0; i < LINES_EXPECTED; i++){
    if (neutron_sep_energy[Z[i]][N[i]] != 0){ 
      sep_energy_plot[Z[i]-1]->SetPoint(n_points[Z[i]-1], N[i], neutron_sep_energy[Z[i]-1][N[i]-1]);
      LD_sep_energy_plot[Z[i]-1]->SetPoint(n_points[Z[i]-1], N[i], neutron_sep_energy_LD[Z[i]-1][N[i]-1]);
      n_points[Z[i]-1] += 1;
    }
  }

  n_equals_z_line->SetLineStyle(7);
  

  TCanvas *drip_canvas  = new TCanvas("drip_canvas", "drip_canvas", 800,600);
  dripline_plots->Add(n_dripline_plot, "l");
  dripline_plots->Add(n_equals_z_line, "l");
  dripline_plots->Draw("a");


  TGraph *binding_energy_plot = new TGraph();
  TGraph *binding_energy_plot_LD = new TGraph();
  TGraph *binding_energy_plot_LD_vol = new TGraph();
  TGraph *binding_energy_plot_LD_volsurf = new TGraph();
  TGraph *binding_energy_plot_LD_volsurfcoul = new TGraph();
  int n_points_binding_energy = 0;
  int n_points_binding_energy_LD = 0;
  int n_points_binding_energy_LD_vol = 0;
  int n_points_binding_energy_LD_volsurf = 0;
  int n_points_binding_energy_LD_volsurfcoul = 0;
  for (int mass = 0; mass < MAX_A; mass++){
    if (binding_energy_per_nucleon[mass] != 0){
      binding_energy_plot->SetPoint(n_points_binding_energy, mass+1, binding_energy_per_nucleon[mass]);
      n_points_binding_energy++;
    }
    if (binding_energy_per_nucleon_LD[mass] != 0){
      binding_energy_plot_LD->SetPoint(n_points_binding_energy_LD, mass+1, binding_energy_per_nucleon_LD[mass]);
      n_points_binding_energy_LD++;
    }

    if (binding_energy_per_nucleon_LD_vol[mass] != 0){
      binding_energy_plot_LD_vol->SetPoint(n_points_binding_energy_LD_vol, mass+1, binding_energy_per_nucleon_LD_vol[mass]);
      n_points_binding_energy_LD_vol++;
    }
    if (binding_energy_per_nucleon_LD_volsurf[mass] != 0){
      binding_energy_plot_LD_volsurf->SetPoint(n_points_binding_energy_LD_volsurf, mass+1, binding_energy_per_nucleon_LD_volsurf[mass]);
      n_points_binding_energy_LD_volsurf++;
    }
    if (binding_energy_per_nucleon_LD_volsurfcoul[mass] != 0){
      binding_energy_plot_LD_volsurfcoul->SetPoint(n_points_binding_energy_LD_volsurfcoul, mass+1, binding_energy_per_nucleon_LD_volsurfcoul[mass]);
      n_points_binding_energy_LD_volsurfcoul++;
    }
  }

  TCanvas *sep_canvas = new TCanvas("sep_canvas", "sep_canvas",800,600);
  sep_energy_plot[Z_OXYGEN-1]->SetMarkerStyle(20);
  sep_energy_plot[Z_CALCIUM-1]->SetMarkerStyle(21);
  sep_energy_plot[Z_LEAD-1]->SetMarkerStyle(22);
  sep_energy_plot[Z_TIN-1]->SetMarkerStyle(34);
  sep_energy_plot[Z_OXYGEN-1]->SetLineColor(kRed);
  sep_energy_plot[Z_CALCIUM-1]->SetLineColor(kBlue);
  sep_energy_plot[Z_LEAD-1]->SetLineColor(kSpring);
  sep_energy_plot[Z_TIN-1]->SetLineColor(kMagenta);

  sep_energy_plot[Z_OXYGEN-1]->SetLineWidth(2);
  sep_energy_plot[Z_CALCIUM-1]->SetLineWidth(2);
  sep_energy_plot[Z_LEAD-1]->SetLineWidth(2);
  sep_energy_plot[Z_TIN-1]->SetLineWidth(2);

  LD_sep_energy_plot[Z_OXYGEN-1]->SetLineColor(kRed);
  LD_sep_energy_plot[Z_CALCIUM-1]->SetLineColor(kBlue);
  LD_sep_energy_plot[Z_LEAD-1]->SetLineColor(kSpring);
  LD_sep_energy_plot[Z_TIN-1]->SetLineColor(kMagenta);

  LD_sep_energy_plot[Z_OXYGEN-1]->SetMarkerStyle(24);
  LD_sep_energy_plot[Z_CALCIUM-1]->SetMarkerStyle(25);
  LD_sep_energy_plot[Z_LEAD-1]->SetMarkerStyle(26);
  LD_sep_energy_plot[Z_TIN-1]->SetMarkerStyle(28);

  LD_sep_energy_plot[Z_OXYGEN-1]->SetLineStyle(7);
  LD_sep_energy_plot[Z_CALCIUM-1]->SetLineStyle(7);
  LD_sep_energy_plot[Z_LEAD-1]->SetLineStyle(7);
  LD_sep_energy_plot[Z_TIN-1]->SetLineStyle(7);

  LD_sep_energy_plot[Z_OXYGEN-1]->SetLineWidth(3);
  LD_sep_energy_plot[Z_CALCIUM-1]->SetLineWidth(3);
  LD_sep_energy_plot[Z_LEAD-1]->SetLineWidth(3);
  LD_sep_energy_plot[Z_TIN-1]->SetLineWidth(3);


  sep_energy_graphs->Add(sep_energy_plot[Z_OXYGEN-1], "lp");
  sep_energy_graphs->Add(sep_energy_plot[Z_CALCIUM-1], "lp");
  sep_energy_graphs->Add(sep_energy_plot[Z_TIN-1], "lp");
  sep_energy_graphs->Add(sep_energy_plot[Z_LEAD-1], "lp");
  sep_energy_graphs->Add(LD_sep_energy_plot[Z_OXYGEN-1], "lp");
  sep_energy_graphs->Add(LD_sep_energy_plot[Z_CALCIUM-1], "lp");
  sep_energy_graphs->Add(LD_sep_energy_plot[Z_TIN-1], "lp");
  sep_energy_graphs->Add(LD_sep_energy_plot[Z_LEAD-1], "lp");

  TLegend *sep_leg = new TLegend (0.8,0.8,0.9,0.9);
  sep_leg->AddEntry(sep_energy_plot[Z_OXYGEN-1], "Oxygen (Z = 8)", "lp");
  sep_leg->AddEntry(sep_energy_plot[Z_CALCIUM-1], "Calcium (Z = 20)", "lp");
  sep_leg->AddEntry(sep_energy_plot[Z_TIN-1], "Tin (Z = 50)", "lp");
  sep_leg->AddEntry(sep_energy_plot[Z_LEAD-1], "Lead (Z = 82)", "lp");
  sep_energy_graphs->Draw("a");
  sep_leg->Draw();

  sep_energy_graphs->GetXaxis()->SetTitle("Neutron number");
  sep_energy_graphs->GetYaxis()->SetTitle("S_{n} (MeV)");

  binding_energy_plot_LD_vol->SetLineColor(kRed);
  binding_energy_plot_LD_volsurf->SetLineColor(kSpring);
  binding_energy_plot_LD_volsurfcoul->SetLineColor(kBlue);

  binding_energy_plot->SetLineWidth(3);
  binding_energy_plot_LD->SetLineWidth(3);
  binding_energy_plot_LD_vol->SetLineWidth(3);
  binding_energy_plot_LD_volsurf->SetLineWidth(3);
  binding_energy_plot_LD_volsurfcoul->SetLineWidth(3);
  TCanvas *be_canvas = new TCanvas("be_canvas", "be_canvas",800,600);
  be_canvas->Divide(1,2);
  be_canvas->cd(1);
  binding_energy_plot->Draw("LA");
  be_canvas->cd(2);
  be_graphs_LD->Add(binding_energy_plot_LD, "l");
  be_graphs_LD->Add(binding_energy_plot_LD_vol, "l");
  be_graphs_LD->Add(binding_energy_plot_LD_volsurf, "l");
  be_graphs_LD->Add(binding_energy_plot_LD_volsurfcoul, "l");
  be_graphs_LD->Draw("A");

  TLegend *be_leg = new TLegend (0.8,0.8,0.9,0.9);
  be_leg->AddEntry(binding_energy_plot_LD, "Total Liquid Drop Binding Energy","l");
  be_leg->AddEntry(binding_energy_plot_LD_vol,"Volume Term", "l");
  be_leg->AddEntry(binding_energy_plot_LD_volsurf,"Volume+Surface Term", "l");
  be_leg->AddEntry(binding_energy_plot_LD_volsurfcoul,"Volume+Surface+Coulomb Term", "l");
  be_leg->Draw();

  binding_energy_plot->GetXaxis()->SetTitle("Mass Number A");
  binding_energy_plot->GetYaxis()->SetTitle("Binding Energy per Nucleon (MeV)");

  be_graphs_LD->GetXaxis()->SetTitle("Mass Number A");
  be_graphs_LD->GetYaxis()->SetTitle("LDM Binding Energy per Nucleon (MeV)");
}






