/*
 * DAnalysis.cpp
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#include "interface/DAnalysis.h"

const double MTOP  = 172.5;
const double MW    = 80.4;
const double Mtlep_mean_0t_  = 175.;
const double Mtlep_sigma_0t_ = 19.;
const double Mthad_mean_0t_  = 177.;
const double Mthad_sigma_0t_ = 16.;

const double Mtlep_mean_1t_  = 175.;
const double Mtlep_sigma_1t_ = 19.;
const double Mthad_mean_1t_  = 173.;
const double Mthad_sigma_1t_ = 15.;

void DAnalysis::analyze(size_t childid /* this info can be used for printouts */){

	/*
	 * This skeleton analyser runs directly on the Delphes output.
	 * It can be used to create histograms directly or a skim.
	 * If a skim is created, a new input configuration will be written automatically
	 * and stored in the output directory together with the ntuples.
	 * The skim can contain delphes objects again or can be flat. This is up
	 * to the user.
	 * Examples for both are given here.
	 *
	 * The same skeleton can be used to read the skim. Please refer to the comments
	 * marked with "==SKIM=="
	 *
	 * These parts are commented, since the code is supposed to work as an example without
	 * modifications on Delphes output directly.
	 */



	/*
	 * Define the branches that are to be considered for the analysis
	 * This branch handler (notice the "d")
	 * is used to run directly in Delphes output.
	 * For skimmed ntuples, see below
	 */
	d_ana::dBranchHandler<Electron> elecs(tree(),"Electron"); //Medium electron
	/*
	 * Other branches might be the following
	 * (for a full list, please inspect the Delphes sample root file with root)
	 * For the Delphes class description, see $DELPHES_PATH/classes/DelphesClasses.h
	 */
	//
	d_ana::dBranchHandler<HepMCEvent>  event(tree(),"Event");
	//d_ana::dBranchHandler<Weight>      weights(tree(),"Weight");
	d_ana::dBranchHandler<GenParticle> genpart(tree(),"Particle");
	//d_ana::dBranchHandler<Jet>         genjet(tree(),"GenJet");
	d_ana::dBranchHandler<Jet>         jet(tree(),"JetPUPPI");
	d_ana::dBranchHandler<Muon>        muontight(tree(),"MuonTight");
	//d_ana::dBranchHandler<Photon>      photon(tree(),"Photon");
	d_ana::dBranchHandler<MissingET>   met(tree(),"PuppiMissingET");
	//d_ana::dBranchHandler<Vertex>      vrtx(tree(),"Vertex");

// 	TString jetCollection = "";
// 	if(getSampleFile().Contains("_0PU")) {jetCollection = "JetAK8";}
// 	else {jetCollection = "JetPUPPIAK8";}
	d_ana::dBranchHandler<Jet> jetAK8(tree(),"JetPUPPIAK8");
	//d_ana::dBranchHandler<Jet> genjetAK8(tree(),"GenJetAK8");
	//std::cout<<"jetCollection="<<jetCollection<<", FName="<<getSampleFile()<<std::endl;

	/* ==SKIM==
	 *
	 * If a skim of the Delphes outout was created in a way indicated
	 * further below, use the tBranchHandler (please notive the "t")
	 * to access vectors of objects...
	 *
	 */
	// d_ana::tBranchHandler<std::vector<Electron> > electrons(tree(),"Electrons");

	/*==SKIM==
	 *
	 * Or an object directly
	 *
	 */
	//d_ana::tBranchHandler<MissingET> met(tree(),"MET");



	/*
	 * Always use this function to add a new histogram (can also be 2D)!
	 * Histograms created this way are automatically added to the output file
	 */
// 	TH1* histo=addPlot(new TH1D("histoname1","histotitle1",100,0,100),"p_{T} [GeV]","N_{e}");
	TFile *YRJECfile = TFile::Open("/afs/cern.ch/work/s/sisagir/private/DAnalysisFW/DAnalysis/DAnalysis_workdir/HL_YR_JEC.root");
	TH1D * htotJESbjets = (TH1D*)YRJECfile->Get("TOTAL_BJES_AntiKt4EMTopo_YR2018");
	TH1D * htotJESljets = (TH1D*)YRJECfile->Get("TOTAL_DIJET_AntiKt4EMTopo_YR2018");

	/*
	 * If (optionally) a skim or a flat ntuple is to be created, please use the following function to initialize
	 * the tree.
	 * The output files will be written automatically, and a config file will be created.
	 */
	TTree* myskim=addTree();
	/*
	 * Add a simple branch to the skim
	 */
	
	//Bool_t doJECup = false;
	//Bool_t doJECdn = false;
	std::cout<<"doJECup="<<doJECup()<<std::endl;
	std::cout<<"doJECdn="<<doJECdn()<<std::endl;
	Bool_t isSignal = getSampleFile().Contains("RSGluonToTTbar_") || getSampleFile().Contains("ZPrimeToTTJets_") || getSampleFile().Contains("CMS_200PU_PhaseII");
	Bool_t isTT  = getSampleFile().Contains("TT_TuneCUETP8M2") || getSampleFile().Contains("TT_Mtt") || getSampleFile().Contains("mg_pp_tt_");
	Bool_t isQCD = getSampleFile().Contains("QCD_Pt-") || getSampleFile().Contains("QCD_Mdijet-") || getSampleFile().Contains("QCD_Flat_") || getSampleFile().Contains("mg_pp_jj_");
	
	Int_t eventNumber=-1;

	Bool_t isSingEl=false;
	Bool_t isSingMu=false;
	
	Double_t metPt =-999;
	Double_t metEta=-999;
	Double_t metPhi=-999;
	Double_t metObjBasedPt =-999;
	Double_t metObjBasedEta=-999;
	Double_t metObjBasedPhi=-999;

	Double_t lepPt =-999;
	Double_t lepEta=-999;
	Double_t lepPhi=-999;
	Double_t lepRelIso=-999;
	Double_t lepAbsIso=-999;
	
	Double_t leadJetPt=-999;
	Double_t leadJetEta=-999;
	Double_t leadJetPhi=-999;
	Double_t leadJetMass=-999;
	Bool_t leadJetBTag=0;
	Double_t subLeadJetPt=-999;
	Double_t subLeadJetEta=-999;
	Double_t subLeadJetPhi=-999;
	Double_t subLeadJetMass=-999;
	Bool_t subLeadJetBTag=0;

// 	std::vector<Double_t> jetPt;
// 	std::vector<Double_t> jetEta;
// 	std::vector<Double_t> jetPhi;
// 	std::vector<Double_t> jetMass;
// 	std::vector<Int_t> jetBTag;
// 	std::vector<Double_t> deltaR_ljets;

	Double_t tlepLeadAK4Pt=-999;
	Double_t tlepLeadAK4Eta=-999;
	Double_t tlepLeadAK4Phi=-999;
	Double_t tlepLeadAK4Mass=-999;
	Bool_t tlepLeadAK4BTag=0;
	Int_t NJetsSel=0;
	
	Double_t minDR_lepJet=999;
	Double_t ptRel_lepJet=-999;

	Double_t WlepPt=-999;
	Double_t WlepEta=-999;
	Double_t WlepPhi=-999;
	Double_t WlepMass=-999;

	Double_t thadPt=-999;
	Double_t thadEta=-999;
	Double_t thadPhi=-999;
	Double_t thadMass=-999;
	Double_t thadChi2=-999;
	Double_t tlepPt=-999;
	Double_t tlepEta=-999;
	Double_t tlepPhi=-999;
	Double_t tlepMass=-999;
	Double_t tlepChi2=-999;

	std::vector<Double_t> jetAK8Pt;
	std::vector<Double_t> jetAK8Eta;
	std::vector<Double_t> jetAK8Phi;
	std::vector<Double_t> jetAK8Mass;
	std::vector<Double_t> jetAK8Tau32;
	std::vector<Double_t> jetAK8SDMass;
	std::vector<Int_t> jetAK8BTag;
	
	Double_t topAK8Pt=-999;
	Double_t topAK8Eta=-999;
	Double_t topAK8Phi=-999;
	Double_t topAK8Mass=-999;
	Double_t topAK8Tau32=-999;
	Double_t topAK8SDMass=-999;
	Bool_t topAK8BTag=0;
	Int_t Ntoptagged=0;
	
	Double_t zpPt=-999;
	Double_t zpEta=-999;
	Double_t zpPhi=-999;
	Double_t zpMass=-999;
	Double_t zpDeltaY=-999;
	Double_t zpDeltaR=-999;
	//Double_t genzpPt=-999;
	//Double_t genzpEta=-999;
	//Double_t genzpPhi=-999;
	Double_t genzpMass=-999;
	TLorentzVector lepton_lv,j0_lv,j1_lv,zp_lv,W_lv,W_lv_r1,W_lv_r2,nu_lv1,nu_lv2,b_lv;

	std::vector<TLorentzVector> genParP4;
	std::vector<Int_t> genParID;
	Double_t genTTorJJMass=-999;
	Double_t genTTorJJPt=-999;
	
	myskim->Branch("eventNumber", &eventNumber);
	myskim->Branch("isSingEl", &isSingEl);
	myskim->Branch("isSingMu", &isSingMu);
	myskim->Branch("lepPt", &lepPt);
	myskim->Branch("lepEta", &lepEta);
	myskim->Branch("lepPhi", &lepPhi);
	myskim->Branch("lepRelIso", &lepRelIso);
	myskim->Branch("lepAbsIso", &lepAbsIso);
	myskim->Branch("metPt", &metPt);
	myskim->Branch("metEta", &metEta);
	myskim->Branch("metPhi", &metPhi);
	myskim->Branch("metObjBasedPt", &metObjBasedPt);
	myskim->Branch("metObjBasedEta", &metObjBasedEta);
	myskim->Branch("metObjBasedPhi", &metObjBasedPhi);
	myskim->Branch("leadJetPt", &leadJetPt);
	myskim->Branch("leadJetEta", &leadJetEta);
	myskim->Branch("leadJetPhi", &leadJetPhi);
	myskim->Branch("leadJetMass", &leadJetMass);
	myskim->Branch("leadJetBTag", &leadJetBTag);
	myskim->Branch("subLeadJetPt", &subLeadJetPt);
	myskim->Branch("subLeadJetEta", &subLeadJetEta);
	myskim->Branch("subLeadJetPhi", &subLeadJetPhi);
	myskim->Branch("subLeadJetMass", &subLeadJetMass);
	myskim->Branch("subLeadJetBTag", &subLeadJetBTag);
	myskim->Branch("tlepLeadAK4Pt", &tlepLeadAK4Pt);
	myskim->Branch("tlepLeadAK4Eta", &tlepLeadAK4Eta);
	myskim->Branch("tlepLeadAK4Phi", &tlepLeadAK4Phi);
	myskim->Branch("tlepLeadAK4Mass", &tlepLeadAK4Mass);
	myskim->Branch("tlepLeadAK4BTag", &tlepLeadAK4BTag);
	myskim->Branch("NJetsSel", &NJetsSel);
	//myskim->Branch("deltaR_ljets", &deltaR_ljets);
	myskim->Branch("minDR_lepJet", &minDR_lepJet);
	myskim->Branch("ptRel_lepJet", &ptRel_lepJet);
	myskim->Branch("WlepPt", &WlepPt);
	myskim->Branch("WlepEta", &WlepEta);
	myskim->Branch("WlepPhi", &WlepPhi);
	myskim->Branch("WlepMass", &WlepMass);
	myskim->Branch("thadPt", &thadPt);
	myskim->Branch("thadEta", &thadEta);
	myskim->Branch("thadPhi", &thadPhi);
	myskim->Branch("thadMass", &thadMass);
	myskim->Branch("thadChi2", &thadChi2);
	myskim->Branch("tlepPt", &tlepPt);
	myskim->Branch("tlepEta", &tlepEta);
	myskim->Branch("tlepPhi", &tlepPhi);
	myskim->Branch("tlepMass", &tlepMass);
	myskim->Branch("tlepChi2", &tlepChi2);
	myskim->Branch("topAK8Pt", &topAK8Pt);
	myskim->Branch("topAK8Eta", &topAK8Eta);
	myskim->Branch("topAK8Phi", &topAK8Phi);
	myskim->Branch("topAK8Mass", &topAK8Mass);
	myskim->Branch("topAK8Tau32", &topAK8Tau32);
	myskim->Branch("topAK8SDMass", &topAK8SDMass);
	myskim->Branch("topAK8BTag", &topAK8BTag);
	myskim->Branch("Ntoptagged", &Ntoptagged);
	myskim->Branch("zpPt", &zpPt);
	myskim->Branch("zpEta", &zpEta);
	myskim->Branch("zpPhi", &zpPhi);
	myskim->Branch("zpMass", &zpMass);
	myskim->Branch("zpDeltaY", &zpDeltaY);
	myskim->Branch("zpDeltaR", &zpDeltaR);
	//myskim->Branch("genzpPt", &genzpPt);
	//myskim->Branch("genzpEta", &genzpEta);
	//myskim->Branch("genzpPhi", &genzpPhi);
	myskim->Branch("genzpMass", &genzpMass);
	myskim->Branch("genTTorJJMass", &genTTorJJMass);
	myskim->Branch("genTTorJJPt", &genTTorJJPt);
	std::vector<TLorentzVector> selectedjets;
	std::vector<Int_t> selectedjetsBTag;
	std::vector<Int_t> selectedjetsFlavor;
	//std::vector<Jet> overlapjets;
	std::vector<TLorentzVector> tlep_jets;
	std::vector<Int_t> tlep_jetsBTag;
	std::vector<Int_t> tlep_jetsFlavor;
	/*
	 * Or store a vector of objects (also possible to store only one object)
	 */
	//std::vector<Electron> skimmedelecs;
	//myskim->Branch("Electrons",&skimmedelecs);



	//std::cout<<"FName="<<getSampleFile()<<" contains _0PU? "<<getSampleFile().Contains("_0PU")<<std::endl;
	Int_t Npass1lep = 0;
	Int_t NpassMET  = 0;
	Int_t Npass2ak4 = 0;
	Int_t Npass2ttag = 0;
	size_t nevents=tree()->entries();
	if(isTestMode())
		nevents/=100;
	for(size_t eventno=0;eventno<nevents;eventno++){
		/*
		 * The following two lines report the status and set the event link
		 * Do not remove!
		 */
		reportStatus(eventno,nevents);
		tree()->setEntry(eventno);



		/*
		 * Begin the event-by-event analysis
		 */
		
		/*
		 * Or to fill the skim
		 */
		//skimmedelecs.clear();
		
		eventNumber = event.at(0)->Number;
		
		Int_t Nels = 0;
		Int_t elIdx = 0;
		for(size_t iel=0;iel<elecs.size();iel++){
			//flat info
			if(elecs.at(iel)->PT < 65) continue;
			if(fabs(elecs.at(iel)->Eta) > 3) continue; //vs. 2.5
			//if(elecs.at(elIdx)->IsolationVarRhoCorr/elecs.at(elIdx)->PT > 2.4) continue;
			Nels++;
			elIdx=iel;
			//or objects
			//skimmedelecs.push_back(*elecs.at(i));
			}

		Int_t Nmus = 0;
		Int_t muIdx = 0;
		for(size_t imu=0;imu<muontight.size();imu++){
			//flat info
			if(muontight.at(imu)->PT < 55) continue;
			if(fabs(muontight.at(imu)->Eta) > 3) continue; //vs. 2.4
			Nmus++;
			muIdx=imu;
			}

		isSingEl=false;
		isSingMu=false;
		lepPt =-999;
		lepEta=-999;
		lepPhi=-999;
		lepRelIso=-999;
		lepAbsIso=-999;
		if(Nels==1 && Nmus==0){
			isSingEl=true;
			lepPt=elecs.at(elIdx)->PT;
			lepEta=elecs.at(elIdx)->Eta;
			lepPhi=elecs.at(elIdx)->Phi;
			lepRelIso=elecs.at(elIdx)->IsolationVarRhoCorr; //Relative isolation
			lepAbsIso=elecs.at(elIdx)->SumPt; //Absolute isolation (SumPt=IsolationVarRhoCorr*pT)
			//lepRelIso=elecs.at(elIdx)->IsolationVarRhoCorr/elecs.at(elIdx)->PT;
			}
		else if(Nels==0 && Nmus==1){
			isSingMu=true;
			lepPt=muontight.at(muIdx)->PT;
			lepEta=muontight.at(muIdx)->Eta;
			lepPhi=muontight.at(muIdx)->Phi;
			lepRelIso=muontight.at(muIdx)->IsolationVarRhoCorr; //Relative isolation
			lepAbsIso=muontight.at(muIdx)->SumPt; //Absolute isolation (SumPt=IsolationVarRhoCorr*pT)
			//lepRelIso=muontight.at(muIdx)->IsolationVarRhoCorr/muontight.at(muIdx)->PT;
			}

		if(!(isSingEl || isSingMu)) continue;
		Npass1lep++;

        float lepM=0;
        if(isSingMu){lepM = 0.105658367;
        }else{lepM = 0.00051099891;}
        lepton_lv.SetPtEtaPhiM(lepPt,lepEta,lepPhi,lepM);
		
		metPt  = met.at(0)->MET;
		metEta = met.at(0)->Eta;
		metPhi = met.at(0)->Phi;

		NJetsSel = 0;
		minDR_lepJet = 1e9; 
		ptRel_lepJet = -999;
		float deltaRljets = 0;
		TLorentzVector correctedJet;
		TLorentzVector correctedMET,objBasedMET;
		selectedjets.clear();
		selectedjetsBTag.clear();
		float metPx=metPt*cos(metPhi);
		float metPy=metPt*sin(metPhi);
		float metObjBasedPx=-lepton_lv.Px();
		float metObjBasedPy=-lepton_lv.Py();
		
		float uncScale = 1;
		int binnum = 0;
		for(size_t ijet=0;ijet<jet.size();ijet++){
			if(jet.at(ijet)->PT < 30) continue;
			if(fabs(jet.at(ijet)->Eta) > 4) continue;
            //Lepton-jet cleaning by subtructing lepton 4v from jet 4v if they aren't separated by DR>0.4
			correctedJet = jet.at(ijet)->P4();
			deltaRljets = lepton_lv.DeltaR(correctedJet);
			if(deltaRljets<0.4) correctedJet = correctedJet-lepton_lv;
			if(correctedJet.Pt() < 30) continue;
			if(fabs(correctedJet.Eta()) > 4) continue;
			if(doJECup() || doJECdn()){
			  	uncScale = 1;
			  	binnum = htotJESbjets->GetXaxis()->FindBin(correctedJet.Pt());
			  	if(correctedJet.Pt()>2999.99){binnum=htotJESbjets->GetSize()-2;}
			  	if(jet.at(ijet)->Flavor==5){uncScale = htotJESbjets->GetBinContent(binnum);}
			    else{uncScale = htotJESljets->GetBinContent(binnum);}
			    if(doJECup()){uncScale = 1 + uncScale;}
			    else if(doJECdn()){uncScale = 1 - uncScale;}
			    correctedJet = correctedJet*uncScale;
			    }
			selectedjets.push_back(correctedJet);
			selectedjetsBTag.push_back(jet.at(ijet)->BTag & (1 << 1)); //Medium WP
			selectedjetsFlavor.push_back(jet.at(ijet)->Flavor);
			deltaRljets = lepton_lv.DeltaR(correctedJet);
			if(deltaRljets<minDR_lepJet){
			  minDR_lepJet = deltaRljets;
			  ptRel_lepJet = lepton_lv.P()*(correctedJet.Vect().Cross(lepton_lv.Vect()).Mag()/correctedJet.P()/lepton_lv.P());
			  }
			metPx += jet.at(ijet)->P4().Px() - correctedJet.Px();
			metPy += jet.at(ijet)->P4().Py() - correctedJet.Py();
			metObjBasedPx=metObjBasedPx-correctedJet.Px();
			metObjBasedPy=metObjBasedPy-correctedJet.Py();
            
            /*//Lepton-jet cleaning by removing jets which aren't separated by DR>0.4 from lepton
 			deltaRljets = lepton_lv.DeltaR(jet.at(ijet)->P4());
			if(deltaRljets<0.4) continue;
			selectedjets.push_back(jet.at(ijet)->P4());
			selectedjetsBTag.push_back(jet.at(ijet)->BTag & (1 << 1)); //Medium WP
			if(deltaRljets<minDR_lepJet){
			  minDR_lepJet = deltaRljets;
			  ptRel_lepJet = lepton_lv.P()*(jet.at(ijet)->P4().Vect().Cross(lepton_lv.Vect()).Mag()/jet.at(ijet)->P4().P()/lepton_lv.P());
			  }*/
			NJetsSel++;
			}

		correctedMET.SetPxPyPzE(metPx,metPy,0,sqrt(metPx*metPx+metPy*metPy));
		metPt  = correctedMET.Pt();
		metEta = correctedMET.Eta();
		metPhi = correctedMET.Phi();
		
		objBasedMET.SetPxPyPzE(metObjBasedPx,metObjBasedPy,0,sqrt(metObjBasedPx*metObjBasedPx+metObjBasedPy*metObjBasedPy));
		metObjBasedPt = objBasedMET.Pt();
		metObjBasedEta = objBasedMET.Eta();
		metObjBasedPhi = objBasedMET.Phi();

		if(metPt<50 && metObjBasedPt<50) continue;
		NpassMET++;
		if(NJetsSel<2) continue;
		if(selectedjets.at(0).Pt()<150) continue;
		if(selectedjets.at(1).Pt()<50) continue;
		Npass2ak4++;

		leadJetPt=selectedjets.at(0).Pt();
		leadJetEta=selectedjets.at(0).Eta();
		leadJetPhi=selectedjets.at(0).Phi();
		leadJetMass=selectedjets.at(0).M();
		leadJetBTag=selectedjetsBTag.at(0);
		subLeadJetPt=selectedjets.at(1).Pt();
		subLeadJetEta=selectedjets.at(1).Eta();
		subLeadJetPhi=selectedjets.at(1).Phi();
		subLeadJetMass=selectedjets.at(1).M();
		subLeadJetBTag=selectedjetsBTag.at(1);
		
		Ntoptagged = 0;
		jetAK8Pt.clear();
		jetAK8Eta.clear();
		jetAK8Phi.clear();
		jetAK8Mass.clear();
		jetAK8Tau32.clear();
		jetAK8SDMass.clear();
		jetAK8BTag.clear();
		TLorentzVector ak8_p4;
		for(size_t ijet=0;ijet<jetAK8.size();ijet++){
			ak8_p4 = jetAK8.at(ijet)->P4();
			if(ak8_p4.Pt() < 400) continue;
			if(fabs(ak8_p4.Eta()) > 4) continue;
			//if(jetAK8.at(ijet)->SoftDroppedJet.M() < 50) continue;
			if(jetAK8.at(ijet)->SoftDroppedJet.M() < 105) continue;
			if(jetAK8.at(ijet)->SoftDroppedJet.M() > 210) continue;
			if(jetAK8.at(ijet)->Tau[2]/jetAK8.at(ijet)->Tau[1] > 0.65) continue;
			if(lepton_lv.DeltaR(ak8_p4)<0.8) continue;
			if(doJECup() || doJECdn()){
			  	uncScale = 1;
			  	binnum = htotJESbjets->GetXaxis()->FindBin(ak8_p4.Pt());
			  	if(ak8_p4.Pt()>2999.99){binnum=htotJESbjets->GetSize()-2;}
			  	if(jetAK8.at(ijet)->Flavor==5){uncScale = htotJESbjets->GetBinContent(binnum);}
			    else{uncScale = htotJESljets->GetBinContent(binnum);}
			    if(doJECup()){uncScale = 1 + uncScale;}
			    else if(doJECdn()){uncScale = 1 - uncScale;}
			    ak8_p4 = ak8_p4*uncScale;
			    }
			jetAK8Pt.push_back(ak8_p4.Pt());
			jetAK8Eta.push_back(ak8_p4.Eta());
			jetAK8Phi.push_back(ak8_p4.Phi());
			jetAK8Mass.push_back(ak8_p4.M());
			jetAK8Tau32.push_back(jetAK8.at(ijet)->Tau[2]/jetAK8.at(ijet)->Tau[1]);
			jetAK8SDMass.push_back(jetAK8.at(ijet)->SoftDroppedJet.M());
			jetAK8BTag.push_back(jetAK8.at(ijet)->BTag);
			Ntoptagged++;
			}
		
		if(Ntoptagged>1) continue; //all-hadronic events
		Npass2ttag++;
					
		tlepLeadAK4Pt=-999;
		tlepLeadAK4Eta=-999;
		tlepLeadAK4Phi=-999;
		tlepLeadAK4Mass=-999;
		tlepLeadAK4BTag=0;
		
		WlepPt=-999;
		WlepEta=-999;
		WlepPhi=-999;
		WlepMass=-999;
		
		thadPt=-999;
		thadEta=-999;
		thadPhi=-999;
		thadMass=-999;
		thadChi2=999;
		
		tlepPt=-999;
		tlepEta=-999;
		tlepPhi=-999;
		tlepMass=-999;
		tlepChi2=999;

		zpPt=-999;
		zpEta=-999;
		zpPhi=-999;
		zpMass=-999;
		zpDeltaY=-999;
		zpDeltaR=-999;
	      
        // ----------------------------------------------------------------------------
        // W --> l nu with mass constraint
        // ----------------------------------------------------------------------------
          	         
        //float metPx = metPt*cos(metPhi);
        //float metPy = metPt*sin(metPhi);
          
        float Dtmp = MW*MW-lepM*lepM+2.0*(lepton_lv.Px()*metPx+lepton_lv.Py()*metPy);
        float Atmp = 4.0*(lepton_lv.Energy()*lepton_lv.Energy()-lepton_lv.Pz()*lepton_lv.Pz());
        float Btmp =-4.0*Dtmp*lepton_lv.Pz();
        float Ctmp = 4.0*lepton_lv.Energy()*lepton_lv.Energy()*metPt*metPt-Dtmp*Dtmp;
        float DETtmp = Btmp*Btmp-4.0*Atmp*Ctmp;
          
        float nuPz_1,nuPz_2;
        if(DETtmp >= 0){ // real roots, two solutions
          nuPz_1 = (-Btmp+TMath::Sqrt(DETtmp))/(2.0*Atmp);
          nuPz_2 = (-Btmp-TMath::Sqrt(DETtmp))/(2.0*Atmp);
        }else{ // complex roots ==> use the real part, only one solution in this case
          nuPz_1 = (-Btmp)/(2.0*Atmp);
          nuPz_2 = (-Btmp)/(2.0*Atmp);
          }
        nu_lv1.SetPxPyPzE(metPx,metPy,nuPz_1,TMath::Sqrt(metPt*metPt+nuPz_1*nuPz_1));
        nu_lv2.SetPxPyPzE(metPx,metPy,nuPz_2,TMath::Sqrt(metPt*metPt+nuPz_2*nuPz_2));

        W_lv_r1 = nu_lv1+lepton_lv;
        W_lv_r2 = nu_lv2+lepton_lv;

	    int bIndex = 0;
	    
	    if(Ntoptagged==0){
	      unsigned int n_neusols = 2;
		  float current_best_disc = std::numeric_limits<float>::infinity();
		  
	      for(unsigned int neu=0; neu < n_neusols; neu++){
			  TLorentzVector wlep_v4;
			  if(neu==0){wlep_v4=W_lv_r1;}
			  else if(neu==1){wlep_v4=W_lv_r2;}
			  
			  unsigned int n_jets = selectedjets.size();
			  if(n_jets>10) n_jets=10; //avoid crashes in events with many jets
			  const unsigned int max_j = pow(3, n_jets);
			  
			  for(unsigned int j=0; j < max_j; j++) {
				 TLorentzVector tophad_v4;
				 tophad_v4.SetPxPyPzE(0,0,0,0);
				 TLorentzVector toplep_v4 = wlep_v4;
				 int nhadjets=0;
				 int nlepjets=0;
				 int num = j;
				 float maxBPt = -1;
				 int bIndexTemp = 0;
				 for(unsigned int ijet=0; ijet < n_jets; ijet++){
					if(num%3==0){
					  tophad_v4 = tophad_v4 + selectedjets.at(ijet);
					  nhadjets++;
					  }
					if(num%3==1){
					  toplep_v4 = toplep_v4 + selectedjets.at(ijet);
					  nlepjets++;
					  if(selectedjets.at(ijet).Pt() > maxBPt){
					  	bIndexTemp=ijet;
					  	maxBPt=selectedjets.at(ijet).Pt();
					  	}
					  }
					num /= 3;
					}
				
				if(nhadjets>0 && nlepjets>0) {
				  const float Mtlep_reco = tophad_v4.M();
				  const float Mthad_reco = toplep_v4.M();
				
				  const double chi2_tlep = pow((Mtlep_reco - Mtlep_mean_0t_) / Mtlep_sigma_0t_, 2);
				  const double chi2_thad = pow((Mthad_reco - Mthad_mean_0t_) / Mthad_sigma_0t_, 2);
				  
				  if(chi2_tlep+chi2_thad < current_best_disc){
				  	j0_lv = tophad_v4;
				  	j1_lv = toplep_v4;
				  	current_best_disc = chi2_tlep+chi2_thad;
				  	thadChi2 = chi2_thad;
				  	tlepChi2 = chi2_tlep;
				  	W_lv = wlep_v4;
				  	bIndex = bIndexTemp;
				  	}
				  }
				} // 3^n_jets jet combinations
			  } // neutrinos
		  }
		
		if(Ntoptagged==1){
	      j0_lv.SetPtEtaPhiM(jetAK8Pt.at(0),jetAK8Eta.at(0),jetAK8Phi.at(0),jetAK8Mass.at(0));
	      TLorentzVector tophad_v4 = j0_lv;

	      tlep_jets.clear();
	      tlep_jetsBTag.clear();
	      for(unsigned int tjet=0; tjet < selectedjets.size(); tjet++){
	      	 if(j0_lv.DeltaR(selectedjets.at(tjet))<1.2) continue;
	      	 tlep_jets.push_back(selectedjets.at(tjet));
	      	 tlep_jetsBTag.push_back(selectedjetsBTag.at(tjet));
	      	 }

		  if(tlep_jets.size()==0){
		  	std::cout<<"There is no AK4 jet with DR>1.2 from t-tagged jet! Skipping ..."<<std::endl;
		  	continue;
		  	}
		  
		  unsigned int n_jets = tlep_jets.size();
		  const unsigned int max_j = pow(2, n_jets);

	      unsigned int n_neusols = 2;
		  float current_best_disc = std::numeric_limits<float>::infinity();
		  
	      for(unsigned int neu=0; neu < n_neusols; neu++){
			  TLorentzVector wlep_v4;
			  if(neu==0){wlep_v4=W_lv_r1;}
			  else if(neu==1){wlep_v4=W_lv_r2;}
			  			  
			  for(unsigned int j=0; j < max_j; j++) {
				 TLorentzVector toplep_v4 = wlep_v4;
				 int nlepjets=0;
				 float maxBPt = -1;
				 int bIndexTemp = 0;
				 for(unsigned int ijet=0; ijet < n_jets; ijet++){
					int jet_topidx = int(j/(pow(2,ijet))) % 2;
					if(jet_topidx == 1){
					  toplep_v4 = toplep_v4 + tlep_jets.at(ijet);
					  nlepjets++;
					  if(tlep_jets.at(ijet).Pt() > maxBPt){
					  	bIndexTemp=ijet;
					  	maxBPt=tlep_jets.at(ijet).Pt();
					  	}
					  }
					}
				
				if(nlepjets>0){
				  const float Mtlep_reco = tophad_v4.M();
				  const float Mthad_reco = toplep_v4.M();
				
				  const double chi2_tlep = pow((Mtlep_reco - Mtlep_mean_1t_) / Mtlep_sigma_1t_, 2);
				  const double chi2_thad = pow((Mthad_reco - Mthad_mean_1t_) / Mthad_sigma_1t_, 2);
				  
				  if(chi2_tlep+chi2_thad < current_best_disc){
				  	j1_lv = toplep_v4;
				  	current_best_disc = chi2_tlep+chi2_thad;
				  	thadChi2 = chi2_thad;
				  	tlepChi2 = chi2_tlep;
				  	W_lv = wlep_v4;
				  	bIndex = bIndexTemp;
				  	}
				  }
				} // 2^n_jets jet combinations
			  } // neutrinos
		  }
		zp_lv = j0_lv+j1_lv;
		
		zpPt = zp_lv.Pt();
		zpEta = zp_lv.Eta();
      	zpPhi = zp_lv.Phi();
      	zpMass = zp_lv.M();
      	zpDeltaY = fabs(j0_lv.Rapidity()-j1_lv.Rapidity());
      	zpDeltaR = j0_lv.DeltaR(j1_lv);

	    if(Ntoptagged==0){
			tlepLeadAK4Pt=selectedjets.at(bIndex).Pt();
			tlepLeadAK4Eta=selectedjets.at(bIndex).Eta();
			tlepLeadAK4Phi=selectedjets.at(bIndex).Phi();
			tlepLeadAK4Mass=selectedjets.at(bIndex).M();
			tlepLeadAK4BTag=selectedjetsBTag.at(bIndex);
	    }else if(Ntoptagged==1){
			tlepLeadAK4Pt=tlep_jets.at(bIndex).Pt();
			tlepLeadAK4Eta=tlep_jets.at(bIndex).Eta();
			tlepLeadAK4Phi=tlep_jets.at(bIndex).Phi();
			tlepLeadAK4Mass=tlep_jets.at(bIndex).M();
			tlepLeadAK4BTag=tlep_jetsBTag.at(bIndex);
			}

      	WlepPt = W_lv.Pt();
      	WlepEta = W_lv.Eta();
      	WlepPhi = W_lv.Phi();
      	WlepMass = W_lv.M();

      	thadPt = j0_lv.Pt();
      	thadEta = j0_lv.Eta();
      	thadPhi = j0_lv.Phi();
      	thadMass = j0_lv.M();

      	tlepPt = j1_lv.Pt();
      	tlepEta = j1_lv.Eta();
      	tlepPhi = j1_lv.Phi();
      	tlepMass = j1_lv.M();

		topAK8Pt=-999;
		topAK8Eta=-999;
		topAK8Phi=-999;
		topAK8Mass=-999;
		topAK8Tau32=-999;
		topAK8SDMass=-999;
		topAK8BTag=0;
		if(Ntoptagged>0){
			topAK8Pt=jetAK8Pt.at(0);
			topAK8Eta=jetAK8Eta.at(0);
			topAK8Phi=jetAK8Phi.at(0);
			topAK8Mass=jetAK8Mass.at(0);
			topAK8Tau32=jetAK8Tau32.at(0);
			topAK8SDMass=jetAK8SDMass.at(0);
			topAK8BTag=jetAK8BTag.at(0);
			}

		//genzpPt=-999;
		//genzpEta=-999;
		//genzpPhi=-999;
		genzpMass=-999;
		if(isSignal){
			for(unsigned int i=0;i<genpart.size();i++){
				if (genpart.at(i)->PID==5100021 || genpart.at(i)->PID==6000047){
					//genzpPt=genpart.at(i)->PT;
					//genzpEta=genpart.at(i)->Eta;
					//genzpPhi=genpart.at(i)->Phi;
					genzpMass=genpart.at(i)->Mass;
					break;
					}
				}
			}
		
		genParP4.clear();
		genParID.clear();
		genTTorJJMass=-999;
		genTTorJJPt=-999;
		if(isTT || isSignal){
			for(unsigned int i=0;i<genpart.size();i++){
				if (fabs(genpart.at(i)->PID)==6 && genpart.at(i)->Status==22){ // "Status 22 : intermediate (intended to have preserved mass)" from http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
					genParP4.push_back(genpart.at(i)->P4());
					genParID.push_back(genpart.at(i)->PID);
					if(genParID.size()>1) break;
					}
				}
			if(genParID.size()!= 2){std::cout<<"WARNING::Found "<<genParID.size()<<" tops!!"<<std::endl;
			}else if(genParID.at(0)*genParID.at(1) > 0){std::cout<<"WARNING::Found 2 tops with the same ID!!";
			}else{
				genTTorJJMass=(genParP4.at(0)+genParP4.at(1)).M();
				genTTorJJPt=(genParP4.at(0)+genParP4.at(1)).Pt();
				}
			}
		
		if(isQCD){
			TLorentzVector genQCDP4;
			genQCDP4.SetPtEtaPhiM(0.,0.,0.,0.);
			for(unsigned int i=0;i<genpart.size();i++){
				if (genpart.at(i)->Status==23){ // "Status 23 : outgoing" from http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
					genQCDP4+=genpart.at(i)->P4();
					genParID.push_back(genpart.at(i)->PID);
					}
				}
			if(genParID.size()== 0){std::cout<<"WARNING::Found no Status 23 particles!!"<<std::endl;
			}else{genTTorJJMass=genQCDP4.M();}
			}

		myskim->Fill();


		/*==SKIM==
		 * Access the branches of the skim
		 */
		//std::vector<Electron> * skimelecs=electrons.content();
		//for(size_t i=0;i<skimelecs->size();i++){
		//	histo->Fill(skimelecs->at(i).PT);
		//}
	}
	
	std::cout<<"Npass1lep /Ntotal    ="<<Npass1lep<<"/"<<nevents<<std::endl;
	std::cout<<"NpassMET  /Npass1lep ="<<NpassMET<<"/"<<Npass1lep<<std::endl;
	std::cout<<"Npass2ak4 /NpassMET  ="<<Npass2ak4<<"/"<<NpassMET<<std::endl;
	std::cout<<"Npass2ttag/Npass2ak4 ="<<Npass2ttag<<"/"<<Npass2ak4<<std::endl;

	/*
	 * Must be called in the end, takes care of thread-safe writeout and
	 * call-back to the parent process
	 */
	processEndFunction();
}



void DAnalysis::postProcess(){
	/*
	 * This function can be used to analyse the output histograms, e.g. extract a signal contribution etc.
	 * The function can also be called directly on an output file with the histograms, if
	 * RunOnOutputOnly = true
	 * is set in the analyser's config file
	 *
	 * This function also represents an example of how the output of the analyser can be
	 * read-back in an external program.
	 * Just include the sampleCollection.h header and follow the procedure below
	 *
	 */

	/*
	 * Here, the input file to the extraction of parameters from the histograms is the output file
	 * of the parallelised analysis.
	 * The sampleCollection class can also be used externally for accessing the output consistently
	 */
// 	d_ana::sampleCollection samplecoll;
// 	samplecoll.readFromFile(getOutPath());
// 
// 	std::vector<TString> alllegends = samplecoll.listAllLegends();

	/*
	 * Example how to process the output.
	 * Usually, one would define the legendname of the histogram to be used here
	 * by hand, e.g. "signal" or "background".
	 * To make this example run in any case, I am using alllegends.at(0), which
	 * could e.g. be the signal legend.
	 *
	 * So in practise, the following would more look like
	 * samplecoll.getHistos("signal");
	 */
// 	if(alllegends.size()>0){
// 		d_ana::histoCollection histos=samplecoll.getHistos(alllegends.at(0));
// 
// 		/*
// 		 * here, the histogram created in the analyze() function is selected and evaluated
// 		 * The histoCollection maintains ownership (you don't need to delete the histogram)
// 		 */
// 		const TH1* myplot=histos.getHisto("histoname1");
// 
// 		std::cout << "(example output): the integral is " << myplot->Integral() <<std::endl;
// 
// 		/*
// 		 * If the histogram is subject to changes, please clone it and take ownership
// 		 */
// 
// 		TH1* myplot2=histos.cloneHisto("histoname1");
// 
// 		/*
// 		 * do something with the histogram
// 		 */
// 
// 		delete myplot2;
// 	}

	/*
	 * do the extraction here.
	 */



}



