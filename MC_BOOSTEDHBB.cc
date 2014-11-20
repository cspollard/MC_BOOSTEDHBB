// -*- C++ -*-
#include "MC_BOOSTEDHBB.hh"

#include <iostream>
#include <map>
#include <string>

#include "Rivet/Tools/Logging.hh"

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/HeavyHadrons.hh"

#include "Rivet/Jet.hh"
#include "Rivet/Projections/FastJets.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/contrib/VariableR.hh"

using std::map;
using std::string;
using namespace Rivet::Cuts;


// global variables. ick
const string ptlab = "$p_T$ / GeV";
const string mlab = "mass / GeV";
const string drlab = "$\\Delta R";
const string etalab = "$\\eta";
const string philab = "$\\phi";
namespace Rivet {

/// @name Analysis methods
//@{

/// Book histograms and initialise projections before the run
void MC_BOOSTEDHBB::init() {
    bookChannel("AllChannels");

    ChargedLeptons clfs(FinalState(-2.5, 2.5, 25*GeV));
    addProjection(clfs, "ChargedLeptons");
    // calo jets constituents
    // TODO
    // don't include high-pt neutrinos or leptons in jets
    // include electrons?
    IdentifiedFinalState nufs(FinalState(-2.5, 2.5, 25*GeV)); //This will look for a range of particle specified by their pid
    nufs.acceptNeutrinos(); //Here we specify what particle we are looking for.In this case any Neutrino.

    MergedFinalState leptonsAndNeutrinos(clfs, nufs);
    addProjection(leptonsAndNeutrinos, "LeptonsAndNeutrinos");

    // track jets constituents
		VetoedFinalState trackParts(ChargedFinalState(-2.5, 2.5, 5*GeV));
		trackParts.addVetoOnThisFinalState(leptonsAndNeutrinos);

    // variable-R jets. With the tracker.
    fastjet::JetDefinition::Plugin *vrPlugTrackRho60 = new fastjet::contrib::VariableRPlugin(60*GeV /* rho < mH */, 0.0, 1.0, fastjet::contrib::VariableRPlugin::AKTLIKE);//Here we create specify the jet definition 
    addProjection(FastJets(trackParts, vrPlugTrackRho60), "AntiKtVRTrackJetsRho60");//Here we add the projection and call it using the string argument. 

    fastjet::JetDefinition::Plugin *vrPlugTrackRho30 = new fastjet::contrib::VariableRPlugin(30*GeV /* rho < mH */, 0.0, 1.0, fastjet::contrib::VariableRPlugin::AKTLIKE);//We do the same again but this time with a Rho=30 
    addProjection(FastJets(trackParts, vrPlugTrackRho30), "AntiKtVRTrackJetsRho30"); 

    fastjet::JetDefinition::Plugin *vrPlugTrackRho120 = new fastjet::contrib::VariableRPlugin(120*GeV /* rho < mH */, 0.0, 1.0, fastjet::contrib::VariableRPlugin::AKTLIKE);//We do the same again but this time with a Rho=120 
    addProjection(FastJets(trackParts, vrPlugTrackRho120), "AntiKtVRTrackJetsRho120"); 

		//This is to look for the b hadrons. We do this to find the exact location of the b-Hadron relative to the centre of the jet.
		addProjection(HeavyHadrons(-2.5,2.5,0.1*GeV), "HeavyHadrons");

		//This part looks for the signature of Higgs event. With missing energy or looking for decay of Vector boson. If not then vetoEvent.
    // TODO
    // minimum pt cutoff?
    MissingMomentum mmfs(FinalState(-4.2, 4.2, 0.5*GeV));
    addProjection(mmfs, "MissingMomentum");

    FinalState fs;
    addProjection(ZFinder(fs, etaIn(-2.5, 2.5) & (pT >= 25*GeV), PID::ELECTRON, 75*GeV, 105*GeV), "ZeeFinder");
    addProjection(ZFinder(fs, etaIn(-2.5, 2.5) & (pT >= 25*GeV), PID::MUON, 75*GeV, 105*GeV), "ZmumuFinder");

    addProjection(WFinder(fs, etaIn(-2.5, 2.5) & (pT > 25*GeV), PID::ELECTRON, 65*GeV, 95*GeV, 25*GeV), "WenuFinder");
    addProjection(WFinder(fs, etaIn(-2.5, 2.5) & (pT > 25*GeV), PID::MUON, 65*GeV, 95*GeV, 25*GeV), "WmunuFinder");

		//Now book histograms to plot
    bookFourMomPair("BHadron-BTrackJet-Rho30");
    bookFourMomPair("BHadron-BTrackJet-Rho60");
    bookFourMomPair("BHadron-BTrackJet-Rho120");

    bookFourMomComp("BHadron-BTrackJet-Rho30");
    bookFourMomComp("BHadron-BTrackJet-Rho60");
    bookFourMomComp("BHadron-BTrackJet-Rho120");

    cutflow = bookHisto1D("cutflow", CUTSLEN, 0, CUTSLEN, "cutflow", "cut", "entries");

    return;
}


/// Perform the per-event analysis
void MC_BOOSTEDHBB::analyze(const Event& event) {
    const double weight = event.weight();

    // reset cut bits.
    for (unsigned int iCut = 0; iCut < CUTSLEN; ++iCut){
			cutBits[iCut] = false;
		}
    cutBits[NONE] = true;

    // find a boson...
    const Particles& zeebosons = applyProjection<ZFinder>(event, "ZeeFinder").bosons();
    const Particles& zmumubosons = applyProjection<ZFinder>(event, "ZmumuFinder").bosons();

    const Particles& wenubosons = applyProjection<WFinder>(event, "WenuFinder").bosons();
    const Particles& wmunubosons =  applyProjection<WFinder>(event, "WmunuFinder").bosons();

    // leptons
    // TODO
    // isolation?
    const Particles& leptons = applyProjection<ChargedLeptons>(event, "ChargedLeptons").particles();

    const Particle& mm = Particle(0, -applyProjection<MissingMomentum>(event, "MissingMomentum").visibleMomentum()); //Note the 4 vector momentum should sum to zero so if visible momentum is none zero is must be balanced in the opposite direction.

    // find vboson
    Particle vboson;
    if (zeebosons.size() && leptons.size() == 2){ //We look for a single Z boson that has decayed into 2 leptons. 
        vboson = zeebosons[0];
        cutBits[ZLL] = true;
    } else if (zmumubosons.size() && leptons.size() == 2) {
        vboson = zmumubosons[0];
        cutBits[ZLL] = true;
    } else if (wenubosons.size() && leptons.size() == 1) {//We look for a single W boson decaying into electron/muon and neutrino. 
        vboson = wenubosons[0];
        cutBits[WLNU] = true;
    } else if (wmunubosons.size() && leptons.size() == 1) {
        vboson = wmunubosons[0];
        cutBits[WLNU] = true;
    } else if (leptons.size() == 0 && mm.pT() > 30*GeV) { //This is looking for a single Z boson decaying into 2 Neutrinos. We look for missing momentum.
        vboson = Particle(23, mm);
        cutBits[ZNUNU] = true;
    }
    // find channel
    string lepchan;
    if (cutBits[ZLL]){
			lepchan = "Zll";
		}else if(cutBits[WLNU]){
			lepchan = "Wlnu";
		}else if (cutBits[ZNUNU]){
			lepchan = "Znunu";
		}else{
			vetoEvent;
		}

    cutBits[VBOSON] = cutBits[ZLL] || cutBits[WLNU] || cutBits[ZNUNU];
		//If we have not vetod then this is a Higgs event and we want to look at the B-Hadrons and B-tagged jets.
    const Jets& antiKtVRTrackJetsRho30 = applyProjection<FastJets>(event, "AntiKtVRTrackJetsRho30").jetsByPt(25*GeV);
		const Jets& antiKtVRTrackJetsBTaggedRho30 = bTagged(antiKtVRTrackJetsRho30);

		const Jets& antiKtVRTrackJetsRho60 = applyProjection<FastJets>(event, "AntiKtVRTrackJetsRho60").jetsByPt(25*GeV);
		const Jets& antiKtVRTrackJetsBTaggedRho60 = bTagged(antiKtVRTrackJetsRho60);

		const Jets& antiKtVRTrackJetsRho120 = applyProjection<FastJets>(event, "AntiKtVRTrackJetsRho120").jetsByPt(25*GeV);
		const Jets& antiKtVRTrackJetsBTaggedRho120 = bTagged(antiKtVRTrackJetsRho120);

		const Particles& bhads = applyProjection<HeavyHadrons>(event, "HeavyHadrons").bHadrons();


		cutBits[ONEBTAGGEDTRACKJETRHO30] = antiKtVRTrackJetsBTaggedRho30.size() == 1; 
		cutBits[TWOBTAGGEDTRACKJETRHO30] = antiKtVRTrackJetsBTaggedRho30.size() == 2; 
		//else this is false and we do not plot.
		cutBits[ONEBTAGGEDTRACKJETRHO60] = antiKtVRTrackJetsBTaggedRho60.size() == 1; 
		cutBits[TWOBTAGGEDTRACKJETRHO60] = antiKtVRTrackJetsBTaggedRho60.size() == 2; 

		cutBits[ONEBTAGGEDTRACKJETRHO120] = antiKtVRTrackJetsBTaggedRho120.size() == 1; 
		cutBits[TWOBTAGGEDTRACKJETRHO120] = antiKtVRTrackJetsBTaggedRho120.size() == 2; 


		//Here we look for b hadrons. We look for b-hadron separate from the jet so we can determine the deltaR between the jet and B-hadron. Note should compare to the nearest jet.  
    if (bhads.size() == 1){
			cutBits[ONEBHADRONSFOUND] = true;
		}
		else if (bhads.size() == 2){
			cutBits[TWOBHADRONSFOUND] = true;
		}else{//We veto the event if no b hadrons are found or if more than 2 are found
			vetoEvent;
		}
    // fill cuts.
    for (int iCut = 0; iCut < CUTSLEN; ++iCut){
			if (cutBits[iCut]){
				cutflow->fill(iCut, weight);
			}
		}
		//Now we want to determine the deltaR between the b-Hadron and the jet which has been b-tagged for a particular channel. It is not important to ensure that this is associated to calo since data should come from Higgs events.
		//We plot first for the small Rho=30 then all others. 
		if(cutBits[ONEBHADRONSFOUND] and cutBits[ONEBTAGGEDTRACKJETRHO30]){
			fillFourMomPair("AllChannels", "BHadron-BTrackJet-Rho30", bhads.at(0).mom(), antiKtVRTrackJetsBTaggedRho30.at(0).mom(), weight);
			fillFourMomComp("AllChannels","BHadron-BTrackJet-Rho30",bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho30.at(0).mom(), weight);
		}
		//This is used to determine the deltaR between the first hadron and it's closest jet. 
		//We also use this to plot if we have a single b-hadron and 2 b tagged jets. We always plot this with the closest jet.
		if((cutBits[TWOBHADRONSFOUND] or cutBits[ONEBHADRONSFOUND]) and cutBits[TWOBTAGGEDTRACKJETRHO30]){
			double dr1 = Rivet::deltaR(bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho30.at(0).mom());
			double dr2 = Rivet::deltaR(bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho30.at(1).mom());
			int trackPosition;
			if(dr1 < dr2){ //So we want to match the first b-hadrons to it's closest jets.
				trackPosition=0;	
			}else{
				trackPosition=1;
			}
			fillFourMomPair("AllChannels", "BHadron-BTrackJet-Rho30", bhads.at(0).mom(), antiKtVRTrackJetsBTaggedRho30.at(trackPosition).mom(), weight);
			fillFourMomComp("AllChannels","BHadron-BTrackJet-Rho30",bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho30.at(trackPosition).mom(), weight);
		}
		//This is used to determine the deltaR between the second hadron and it's closest jet. 
		if(cutBits[TWOBHADRONSFOUND] and cutBits[TWOBTAGGEDTRACKJETRHO30]){
			double dr1 = Rivet::deltaR(bhads.at(1).mom(),antiKtVRTrackJetsBTaggedRho30.at(0).mom());
			double dr2 = Rivet::deltaR(bhads.at(1).mom(),antiKtVRTrackJetsBTaggedRho30.at(1).mom());
			int trackPosition;
			if(dr1 < dr2){ //So we want to match the first b-hadrons to it's closest jets.
				trackPosition=0;	
			}else{
				trackPosition=1;
			}
			fillFourMomPair("AllChannels", "BHadron-BTrackJet-Rho30", bhads.at(1).mom(), antiKtVRTrackJetsBTaggedRho30.at(trackPosition).mom(), weight);
			fillFourMomComp("AllChannels","BHadron-BTrackJet-Rho30",bhads.at(1).mom(),antiKtVRTrackJetsBTaggedRho30.at(trackPosition).mom(), weight);
		}

		if(cutBits[ONEBHADRONSFOUND] and cutBits[ONEBTAGGEDTRACKJETRHO60]){
			fillFourMomPair("AllChannels", "BHadron-BTrackJet-Rho60", bhads.at(0).mom(), antiKtVRTrackJetsBTaggedRho60.at(0).mom(), weight);
			fillFourMomComp("AllChannels","BHadron-BTrackJet-Rho60",bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho60.at(0).mom(), weight);

		}
		//This is used to determine the deltaR between the first hadron and it's closest jet. 
		if((cutBits[TWOBHADRONSFOUND] or cutBits[ONEBHADRONSFOUND]) and cutBits[TWOBTAGGEDTRACKJETRHO60]){
			double dr1 = Rivet::deltaR(bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho60.at(0).mom());
			double dr2 = Rivet::deltaR(bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho60.at(1).mom());
			int trackPosition;
			if(dr1 < dr2){ //So we want to match the first b-hadrons to it's closest jets.
				trackPosition=0;	
			}else{
				trackPosition=1;
			}
			fillFourMomPair("AllChannels", "BHadron-BTrackJet-Rho60", bhads.at(0).mom(), antiKtVRTrackJetsBTaggedRho60.at(trackPosition).mom(), weight);
			fillFourMomComp("AllChannels","BHadron-BTrackJet-Rho60",bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho60.at(trackPosition).mom(), weight);

		}
		//This is used to determine the deltaR between the second hadron and it's closest jet. 
		if(cutBits[TWOBHADRONSFOUND] and cutBits[TWOBTAGGEDTRACKJETRHO60]){
			double dr1 = Rivet::deltaR(bhads.at(1).mom(),antiKtVRTrackJetsBTaggedRho60.at(0).mom());
			double dr2 = Rivet::deltaR(bhads.at(1).mom(),antiKtVRTrackJetsBTaggedRho60.at(1).mom());
			int trackPosition;
			if(dr1 < dr2){ //So we want to match the first b-hadrons to it's closest jets.
				trackPosition=0;	
			}else{
				trackPosition=1;
			}
			fillFourMomPair("AllChannels", "BHadron-BTrackJet-Rho60", bhads.at(1).mom(), antiKtVRTrackJetsBTaggedRho60.at(trackPosition).mom(), weight);
			fillFourMomComp("AllChannels","BHadron-BTrackJet-Rho60",bhads.at(1).mom(),antiKtVRTrackJetsBTaggedRho60.at(trackPosition).mom(), weight);
		}

		if(cutBits[ONEBHADRONSFOUND] and cutBits[ONEBTAGGEDTRACKJETRHO120]){
			fillFourMomPair("AllChannels", "BHadron-BTrackJet-Rho120", bhads.at(0).mom(), antiKtVRTrackJetsBTaggedRho120.at(0).mom(), weight);
			fillFourMomComp("AllChannels","BHadron-BTrackJet-Rho120",bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho120.at(0).mom(), weight);

		}
		//This is used to determine the deltaR between the first hadron and it's closest jet. 
		if((cutBits[TWOBHADRONSFOUND] or cutBits[ONEBHADRONSFOUND]) and cutBits[TWOBTAGGEDTRACKJETRHO120]){
			double dr1 = Rivet::deltaR(bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho120.at(0).mom());
			double dr2 = Rivet::deltaR(bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho120.at(1).mom());
			int trackPosition;
			if(dr1 < dr2){ //So we want to match the first b-hadrons to it's closest jets.
				trackPosition=0;	
			}else{
				trackPosition=1;
			}
			fillFourMomPair("AllChannels", "BHadron-BTrackJet-Rho120", bhads.at(0).mom(), antiKtVRTrackJetsBTaggedRho120.at(trackPosition).mom(), weight);
			fillFourMomComp("AllChannels","BHadron-BTrackJet-Rho120",bhads.at(0).mom(),antiKtVRTrackJetsBTaggedRho120.at(trackPosition).mom(), weight);
		}
		//This is used to determine the deltaR between the second hadron and it's closest jet. 
		if(cutBits[TWOBHADRONSFOUND] and cutBits[TWOBTAGGEDTRACKJETRHO120]){
			double dr1 = Rivet::deltaR(bhads.at(1).mom(),antiKtVRTrackJetsBTaggedRho120.at(0).mom());
			double dr2 = Rivet::deltaR(bhads.at(1).mom(),antiKtVRTrackJetsBTaggedRho120.at(1).mom());
			int trackPosition;
			if(dr1 < dr2){ //So we want to match the first b-hadrons to it's closest jets.
				trackPosition=0;	
			}else{
				trackPosition=1;
			}
			fillFourMomPair("AllChannels", "BHadron-BTrackJet-Rho120", bhads.at(1).mom(), antiKtVRTrackJetsBTaggedRho120.at(trackPosition).mom(), weight);
			fillFourMomComp("AllChannels","BHadron-BTrackJet-Rho120",bhads.at(1).mom(),antiKtVRTrackJetsBTaggedRho120.at(trackPosition).mom(), weight);
		}
    return;
}


/// Normalise histograms etc., after the run
void MC_BOOSTEDHBB::finalize() {

    // normalize to 1/fb
    double norm = 1000*crossSection()/sumOfWeights();
    for (map< string, map<string, map<string, Histo1DPtr> > >::iterator p = histos1D.begin(); p != histos1D.end(); ++p) {
        for (map<string, map<string, Histo1DPtr> >::iterator q = p->second.begin(); q != p->second.end(); ++q) {
            for (map<string, Histo1DPtr>::iterator r = q->second.begin(); r != q->second.end(); ++r) {
                r->second->scaleW(norm); // norm to cross section
            }
        }
    }


    for (map< string, map<string, map<string, Histo2DPtr> > >::iterator p = histos2D.begin(); p != histos2D.end(); ++p) {
        for (map<string, map<string, Histo2DPtr> >::iterator q = p->second.begin(); q != p->second.end(); ++q) {
            for (map<string, Histo2DPtr>::iterator r = q->second.begin(); r != q->second.end(); ++r) {
                r->second->scaleW(norm); // norm to cross section
            }
        }
    }


    cutflow->scaleW(norm);


    return;
}


void MC_BOOSTEDHBB::bookChannel(const string& channel) {

    channels.push_back(channel);
    histos1D[channel] = map<string, map<string, Histo1DPtr> >();
    histos2D[channel] = map<string, map<string, Histo2DPtr> >();

    return;
}

//Book 1D histogram
Histo1DPtr MC_BOOSTEDHBB::bookHisto(const string& name, const string& title,
        const string& xlabel, int nxbins, double xmin, double xmax) {

    double xbinwidth = (xmax - xmin)/nxbins;

    char buff[100];
    sprintf(buff, "events / %.2f", xbinwidth);
    string ylabel = buff;

    return bookHisto1D(name, nxbins, xmin, xmax, title, xlabel, ylabel);
}


//Book 2D histogram.
Histo2DPtr MC_BOOSTEDHBB::bookHisto(const string& name, const string& title,
        const string& xlabel, int nxbins, double xmin, double xmax,
        const string& ylabel, int nybins, double ymin, double ymax) {

    double xbinwidth = (xmax - xmin)/nxbins;
    double ybinwidth = (ymax - ymin)/nybins;

    char buff[100];
    sprintf(buff, "events / %.2f / %.2f", xbinwidth, ybinwidth);
    string zlabel = buff;

    return bookHisto2D(name, nxbins, xmin, xmax, nybins, ymin, ymax, title, xlabel, ylabel, zlabel);
}
//Book 1D profile here.
Profile1DPtr MC_BOOSTEDHBB::bookProfile(const string& name, const string& title, const string& xlabel, int nxbins, double xmin, double xmax){

	double xbinwidth = (xmax - xmin)/nxbins;

	char buff[100];
	sprintf(buff, "events / %.2f", xbinwidth);
	string ylabel = buff;

	return bookProfile1D(name, nxbins, xmin, xmax, title, xlabel, ylabel);
}

//Here book variables which are usually plotted for all projections.
void MC_BOOSTEDHBB::bookFourMom(const string& name) {
    MSG_DEBUG("Booking " << name << " histograms.");

    foreach (const string& chan, channels) {
        histos1D[chan][name]["pt"] = bookHisto(chan + "_" + name + "_pt", name, ptlab, 25, 0, 2000*GeV);
        histos1D[chan][name]["eta"] = bookHisto(chan + "_" + name + "_eta", name, "$\\eta$", 25, -5, 5);
        histos1D[chan][name]["m"] = bookHisto(chan + "_" + name + "_m", name, mlab, 25, 0, 1000*GeV);

        histos2D[chan][name]["m_vs_pt"] = bookHisto(chan + "_" + name + "_m_vs_pt", name,
                ptlab, 25, 0, 2000*GeV,
                mlab, 25, 0, 1000*GeV);
    }

    return;
}


void MC_BOOSTEDHBB::bookFourMomPair(const string& name) {
    // pairs of particles also are "particles"
//    bookFourMom(name); Do not need to plot all this information.

    foreach (const string& chan, channels) {
        histos2D[chan][name]["dr_vs_ptBHad"] = bookHisto(chan + "_" + name + "_dr_vs_ptBHad", name,
                ptlab, 25, 0, 600*GeV,
                drlab, 100, 0, 0.05);
        histos2D[chan][name]["dr_vs_ptTrackJet"] = bookHisto(chan + "_" + name + "_dr_vs_ptTrackJet", name,
                ptlab, 25, 0, 600*GeV,
                drlab, 100, 0, 0.05);
        profile1D[chan][name]["dr_vs_ptBHad"] = bookProfile(chan + "_" + name + "_dr_vs_ptBHad", name,
                ptlab, 25, 0, 600*GeV);
        profile1D[chan][name]["dr_vs_ptTrackJet"] = bookProfile(chan + "_" + name + "_dr_vs_ptTrackJet", name,
                ptlab, 25, 0, 600*GeV);
 
    }

    return;
}

//Book a lot of useful comparison histograms.
void MC_BOOSTEDHBB::bookFourMomComp(const string& name) {

	foreach (const string& chan, channels) {
		histos1D[chan][name]["dr"] = bookHisto(chan + "_" + name + "_dr", drlab, name, 25, 0, 0.5);
	}

    return;
}


void MC_BOOSTEDHBB::bookFourMomColl(const string& name) {
    bookFourMom(name);

     bookFourMom(name + "0");
     bookFourMom(name + "1");

    foreach (const string& chan, channels)
        histos1D[chan][name]["n"] = bookHisto(chan + "_" + name + "_n", "multiplicity", "", 10, 0, 10);

    return;
}


void MC_BOOSTEDHBB::fillFourMom(const string& channel, const string& name, const FourMomentum& p, double weight) {
    MSG_DEBUG("Filling " << name << " histograms");

    histos1D[channel][name]["pt"]->fill(p.pT(), weight);
    histos1D[channel][name]["eta"]->fill(p.eta(), weight);
    histos1D[channel][name]["m"]->fill(p.mass(), weight);
    histos2D[channel][name]["m_vs_pt"]->fill(p.pT(), p.mass(), weight);

    return;
}


void MC_BOOSTEDHBB::fillFourMomPair(const string& channel, const string& name, const FourMomentum& p1, const FourMomentum& p2, double weight) {
//    fillFourMom(channel, name, p1 + p2, weight);

    double dr = Rivet::deltaR(p1, p2);

    histos2D[channel][name]["dr_vs_ptBHad"]->fill(p1.pT(), dr);
    histos2D[channel][name]["dr_vs_ptTrackJet"]->fill(p2.pT(), dr);
    profile1D[channel][name]["dr_vs_ptBHad"]->fill(p1.pT(), dr);
		profile1D[channel][name]["dr_vs_ptTrackJet"]->fill(p2.pT(), dr);
    return;
}


void MC_BOOSTEDHBB::fillFourMomComp(const string& channel, const string& name, const FourMomentum& p1, const FourMomentum& p2, double weight) {

    double dr = Rivet::deltaR(p1, p2);

    histos1D[channel][name]["dr"]->fill(dr, weight);

    return;
}


template <class T>
void MC_BOOSTEDHBB::fillFourMomColl(const string& channel, const string& name, const vector<T>& ps, double weight) {

    MSG_DEBUG("Filling " << ps.size() << " members of collection " << name);
    histos1D[channel][name]["n"]->fill(ps.size(), weight);

		foreach (const T& p, ps)
        fillFourMom(channel, name, p.mom(), weight);

    return;
}


Jets MC_BOOSTEDHBB::bTagged(const Jets& js) {
    Jets bjs;
    foreach(const Jet& j, js)
        // TODO
        // update for new rivet version
        if (j.bTags().size()) bjs.push_back(j);

    return bjs;
}

//@}

} // Rivet
