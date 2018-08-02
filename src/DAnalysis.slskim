/*
 * DAnalysis.cpp
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#include "interface/DAnalysis.h"

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
	//d_ana::dBranchHandler<Weight>      weights(tree(),"Weight");d
	//d_ana::dBranchHandler<GenParticle> genpart(tree(),"Particle");
	//d_ana::dBranchHandler<Jet>         genjet(tree(),"GenJet");
	d_ana::dBranchHandler<Jet>         jet(tree(),"Jet");
	d_ana::dBranchHandler<Muon>        muontight(tree(),"MuonTight");
	//d_ana::dBranchHandler<Photon>      photon(tree(),"Photon");
	d_ana::dBranchHandler<MissingET>   met(tree(),"MissingET");
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


	/*
	 * If (optionally) a skim or a flat ntuple is to be created, please use the following function to initialize
	 * the tree.
	 * The output files will be written automatically, and a config file will be created.
	 */
	TTree* myskim=addTree();
	/*
	 * Add a simple branch to the skim
	 */
	
	Bool_t isSingEl=false;
	Bool_t isSingMu=false;
	Bool_t isAllHad=false;
	
	Double_t lepPt =-999;
	Double_t lepEta=-999;
	Double_t lepPhi=-999;

	Double_t metPt =-999;
	Double_t metEta=-999;
	Double_t metPhi=-999;

	std::vector<Double_t> jetPt;
	std::vector<Double_t> jetEta;
	std::vector<Double_t> jetPhi;
	std::vector<Double_t> jetMass;
	std::vector<Double_t> jetBTag;
	
	Double_t jetAK8Pt=-999;
	Double_t jetAK8Eta=-999;
	Double_t jetAK8Phi=-999;
	Double_t jetAK8Mass=-999;
	Double_t jetAK8Tau1=-999;
	Double_t jetAK8Tau2=-999;
	Double_t jetAK8Tau3=-999;
	Double_t jetAK8SDMass=-999;
	Double_t jetAK8BTag=-999;
	
	myskim->Branch("isSingEl", &isSingEl);
	myskim->Branch("isSingMu", &isSingMu);
	myskim->Branch("isAllHad", &isAllHad);
	myskim->Branch("lepPt", &lepPt);
	myskim->Branch("lepEta", &lepEta);
	myskim->Branch("lepPhi", &lepPhi);
	myskim->Branch("metPt", &metPt);
	myskim->Branch("metEta", &metEta);
	myskim->Branch("metPhi", &metPhi);
	myskim->Branch("jetPt", &jetPt);
	myskim->Branch("jetEta", &jetEta);
	myskim->Branch("jetPhi", &jetPhi);
	myskim->Branch("jetMass", &jetMass);
	myskim->Branch("jetBTag", &jetBTag);
	myskim->Branch("jetAK8Pt", &jetAK8Pt);
	myskim->Branch("jetAK8Eta", &jetAK8Eta);
	myskim->Branch("jetAK8Phi", &jetAK8Phi);
	myskim->Branch("jetAK8Mass", &jetAK8Mass);
	myskim->Branch("jetAK8Tau1", &jetAK8Tau1);
	myskim->Branch("jetAK8Tau2", &jetAK8Tau2);
	myskim->Branch("jetAK8Tau3", &jetAK8Tau3);
	myskim->Branch("jetAK8SDMass", &jetAK8SDMass);
	myskim->Branch("jetAK8BTag", &jetAK8BTag);
	/*
	 * Or store a vector of objects (also possible to store only one object)
	 */
	//std::vector<Electron> skimmedelecs;
	//myskim->Branch("Electrons",&skimmedelecs);



	//std::cout<<"FName="<<getSampleFile()<<" contains _0PU? "<<getSampleFile().Contains("_0PU")<<std::endl;
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
		 
		isSingEl=false;
		isSingMu=false;
		isAllHad=false;
		lepPt =-999;
		lepEta=-999;
		lepPhi=-999;
		metPt =-999;
		metEta=-999;
		metPhi=-999;
		jetAK8Pt=-999;
		jetAK8Eta=-999;
		jetAK8Phi=-999;
		jetAK8Mass=-999;
		jetAK8Tau1=-999;
		jetAK8Tau2=-999;
		jetAK8Tau3=-999;
		jetAK8SDMass=-999;
		jetAK8BTag=-999;
		
		/*
		 * Or to fill the skim
		 */
		//skimmedelecs.clear();
		Int_t Nels = 0;
		Int_t elIdx = 0;
		for(size_t iel=0;iel<elecs.size();iel++){
			//flat info
			if(elecs.at(iel)->PT < 20) continue;
			if(fabs(elecs.at(iel)->Eta) > 2.4) continue;
			Nels++;
			elIdx=iel;
			//or objects
			//skimmedelecs.push_back(*elecs.at(i));
			}

		Int_t Nmus = 0;
		Int_t muIdx = 0;
		for(size_t imu=0;imu<muontight.size();imu++){
			//flat info
			if(muontight.at(imu)->PT < 20) continue;
			if(fabs(muontight.at(imu)->Eta) > 2.4) continue;
			Nmus++;
			muIdx=imu;
			}

		if(Nels==0 && Nmus==0){isAllHad=true;} 
		else if(Nels==1 && Nmus==0){
			isSingEl=true;
			lepPt=elecs.at(elIdx)->PT;
			lepEta=elecs.at(elIdx)->Eta;
			lepPhi=elecs.at(elIdx)->Phi;
			}
		else if(Nels==0 && Nmus==1){
			isSingMu=true;
			lepPt=muontight.at(muIdx)->PT;
			lepEta=muontight.at(muIdx)->Eta;
			lepPhi=muontight.at(muIdx)->Phi;
			}
		
		if(!(isSingEl || isSingMu)) continue;
		if(met.at(0)->MET<30) continue;
		
		metPt  = met.at(0)->MET;
		metEta = met.at(0)->Eta;
		metPhi = met.at(0)->Phi;

		Int_t Njets = 0;
		jetPt.clear();
		jetEta.clear();
		jetPhi.clear();
		jetMass.clear();
		jetBTag.clear();
		for(size_t ijet=0;ijet<jet.size();ijet++){
			if(jet.at(ijet)->PT < 30) continue;
			if(fabs(jet.at(ijet)->Eta) > 4) continue;
			jetPt.push_back(jet.at(ijet)->PT);
			jetEta.push_back(jet.at(ijet)->Eta);
			jetPhi.push_back(jet.at(ijet)->Phi);
			jetMass.push_back(jet.at(ijet)->Mass);
			jetBTag.push_back(jet.at(ijet)->BTag);
			Njets++;
			}
		
		if(Njets==0) continue;
		
		Int_t NAK8jets = 0;
		Int_t ak8Index = 0;
		for(size_t ijet=0;ijet<jetAK8.size();ijet++){
			if(jetAK8.at(ijet)->PT < 400) continue;
			if(fabs(jetAK8.at(ijet)->Eta) > 4) continue;
			if(jetAK8.at(ijet)->SoftDroppedJet.M() < 50) continue;
			NAK8jets++;
			ak8Index = ijet;
			}
		if(NAK8jets==0) continue;
		
		jetAK8Pt=jetAK8.at(ak8Index)->PT;
		jetAK8Eta=jetAK8.at(ak8Index)->Eta;
		jetAK8Phi=jetAK8.at(ak8Index)->Phi;
		jetAK8Mass=jetAK8.at(ak8Index)->Mass;
		jetAK8Tau1=jetAK8.at(ak8Index)->Tau[0];
		jetAK8Tau2=jetAK8.at(ak8Index)->Tau[1];
		jetAK8Tau3=jetAK8.at(ak8Index)->Tau[2];
		jetAK8SDMass=jetAK8.at(ak8Index)->SoftDroppedJet.M();
		jetAK8BTag=jetAK8.at(ak8Index)->BTag;

		myskim->Fill();


		/*==SKIM==
		 * Access the branches of the skim
		 */
		//std::vector<Electron> * skimelecs=electrons.content();
		//for(size_t i=0;i<skimelecs->size();i++){
		//	histo->Fill(skimelecs->at(i).PT);
		//}
	}


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



