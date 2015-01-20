// -*- C++ -*-
#include "MC_BOOSTEDHBB.hh"

#include <iostream>
#include <map>
#include <string>

#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"

#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/FastJets.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/contrib/VariableR.hh"

using std::map;
using std::string;
using namespace Rivet::Cuts;


// global variables. ick
const string ptlab = "$p_T$ / GeV";
const string mlab = "mass / GeV";
const string drlab = "$\\Delta R$";
const string etalab = "$\\eta$";
const string philab = "$\\phi$";



namespace Rivet {



    fastjet::JetDefinition::Plugin *aktVRPlugin(double rho, double rmin, double rmax) {
        return new fastjet::contrib::VariableRPlugin(rho, rmin, rmax,
                fastjet::contrib::VariableRPlugin::AKTLIKE);
    }


    int drMinJetIdx(const Particle& p, const Jets& jets) {

        double dr;
        double drmin = -1;
        int drminidx = -1;
        for (unsigned int iJet = 0; iJet < jets.size(); iJet++) {
            dr = Rivet::deltaR(p.mom(), jets[iJet].mom());
            if (dr < drmin || drminidx < 0) {
                drmin = dr;
                drminidx = iJet;
            }
        }

        return drminidx;
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void MC_BOOSTEDHBB::init() {

        collections.clear();

        bookChannel("Rho10Min00Max1");
        bookChannel("Rho20Min00Max1");
        bookChannel("Rho30Min00Max1");

        bookChannel("AKTTrack02");
        bookChannel("AKTTrack03");
        bookChannel("AKTCalo04");


        // prepare the jet collections with their minimum pt cuts.
        collections.push_back(make_pair("Rho10Min00Max1", 10*GeV));
        collections.push_back(make_pair("Rho20Min00Max1", 10*GeV));
        collections.push_back(make_pair("Rho30Min00Max1", 10*GeV));
        collections.push_back(make_pair("AKTTrack02", 10*GeV));
        collections.push_back(make_pair("AKTTrack03", 10*GeV));
        collections.push_back(make_pair("AKTCalo04", 25*GeV));


        // calo jets constituents
        FinalState caloParts(-2.5, 2.5, 0.5*GeV);

        // track jets constituents
        ChargedFinalState trackParts(-2.5, 2.5, 0.5*GeV);


        // variable-R jet projections
        addProjection(FastJets(trackParts, aktVRPlugin(10*GeV, 0, 1)),
                "Rho10Min00Max1");
        addProjection(FastJets(trackParts, aktVRPlugin(20*GeV, 0, 1)),
                "Rho20Min00Max1");
        addProjection(FastJets(trackParts, aktVRPlugin(30*GeV, 0, 1)),
                "Rho30Min00Max1");


        // conventional jet projections
        addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.2), "AKTTrack02");
        addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.3), "AKTTrack03");
        addProjection(FastJets(caloParts, FastJets::ANTIKT, 0.4), "AKTCalo04");


        // projection to find b-hadrons
        addProjection(HeavyHadrons(-2.5, 2.5, 5*GeV), "HeavyHadrons");

        //Now book histograms to plot
        bookFourMom("GABHad");
        bookFourMom("GABHadJet");
        bookFourMomComp("GABHad", "Jet");
        bookFourMomAllChannels("AllBHad");
        return;
    }


    /// Perform the per-event analysis
    void MC_BOOSTEDHBB::analyze(const Event& event) {
        const Particles& bhads = applyProjection<HeavyHadrons>(event, "HeavyHadrons").bHadrons();

        if (!bhads.size()) {
            vetoEvent;
        }

        foreach (const Particle bHad, bhads)
            fillFourMomAllChannels("AllBHad", bHad, event.weight());

        fillGABHadHists(event, collections);

        return;
    }


    /// Normalise histograms etc., after the run
    void MC_BOOSTEDHBB::finalize() {

        // normalize to 1/fb
        // We must remember to normalise the graphs to the total number of events measured in fb. 
        // Note the weights are used since we produce some very unlikely events more often than predicted so we do not have to run for a very long time.
        double norm = 1000*crossSection()/sumOfWeights();
        for (map< string, map<string, map<string, Histo1DPtr> > >::iterator p = histos1D.begin(); p != histos1D.end(); ++p) {
            for (map<string, map<string, Histo1DPtr> >::iterator q = p->second.begin(); q != p->second.end(); ++q) {
                for (map<string, Histo1DPtr>::iterator r = q->second.begin(); r != q->second.end(); ++r) {
                    r->second->scaleW(norm); // norm to cross section
                }
            }
        }


        for (map< string, map<string, map<string, Profile1DPtr> > >::iterator p = profiles1D.begin(); p != profiles1D.end(); ++p) {
            for (map<string, map<string, Profile1DPtr> >::iterator q = p->second.begin(); q != p->second.end(); ++q) {
                for (map<string, Profile1DPtr>::iterator r = q->second.begin(); r != q->second.end(); ++r) {
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
        //Here we normalise all B-Had graph with pT.
        histos1DAllChannels["AllBHad"]["pt"]->scaleW(norm); 
        histos1DAllChannels["AllBHad"]["eta"]->scaleW(norm);


        return;
    }


    void MC_BOOSTEDHBB::bookChannel(const string& channel) {

        channels.push_back(channel);
        histos1D[channel] = map<string, map<string, Histo1DPtr> >();
        histos2D[channel] = map<string, map<string, Histo2DPtr> >();

        return;
    }


    Histo1DPtr MC_BOOSTEDHBB::bookHisto(const string& label, const string& title,
            const string& xlabel, int nxbins, double xmin, double xmax) {

        double xbinwidth = (xmax - xmin)/nxbins;

        char buff[100];
        sprintf(buff, "events / %.2f", xbinwidth);
        string ylabel = buff;

        return bookHisto1D(label, nxbins, xmin, xmax, title, xlabel, ylabel);
    }


    Histo2DPtr MC_BOOSTEDHBB::bookHisto(const string& label, const string& title,
            const string& xlabel, int nxbins, double xmin, double xmax,
            const string& ylabel, int nybins, double ymin, double ymax) {

        double xbinwidth = (xmax - xmin)/nxbins;
        double ybinwidth = (ymax - ymin)/nybins;

        char buff[100];
        sprintf(buff, "events / %.2f / %.2f", xbinwidth, ybinwidth);
        string zlabel = buff;

        return bookHisto2D(label, nxbins, xmin, xmax, nybins, ymin, ymax, title, xlabel, ylabel, zlabel);
    }


    Profile1DPtr MC_BOOSTEDHBB::bookProfile(const string& label, const string& title,
            const string& xlabel, int nxbins, double xmin, double xmax,
            const string& ylabel) {

        return bookProfile1D(label, nxbins, xmin, xmax, title, xlabel, ylabel);
    }


    void MC_BOOSTEDHBB::bookFourMom(const string& label) {
        MSG_DEBUG("Booking " << label << " histograms.");

        foreach (const string& chan, channels) {
            histos1D[chan][label]["pt"] = bookHisto(chan + "_" + label + "_pt", label, ptlab, 25, 0, 2000*GeV);
            histos1D[chan][label]["eta"] = bookHisto(chan + "_" + label + "_eta", label, "$\\eta$", 25, -5, 5);
        }

        return;
    }

    void MC_BOOSTEDHBB::bookFourMomAllChannels(const string& label) {
        MSG_DEBUG("Booking " << label << " histograms.");

        histos1DAllChannels[label]["pt"] = bookHisto(label + "_pt", label, ptlab, 25, 0, 2000*GeV);
        histos1DAllChannels[label]["eta"] = bookHisto(label + "_eta", label, "$\\eta$", 25, -5, 5);

        return;
    }

    void MC_BOOSTEDHBB::bookFourMomPair(const string& label1,
            const string& label2) {

        const string label = label1 + "_" + label2;

        // pairs of particles also are "particles"
        bookFourMom(label);

        foreach (const string& chan, channels) {
            // extra histograms for pairs of particles
            histos1D[chan][label]["dr"] = bookHisto(chan + "_" + label + "_dr", drlab, "", 25, 0, 5);

            histos2D[chan][label]["dr_vs_pt"] = bookHisto(chan + "_" + label + "_dr_vs_pt", label,
                    ptlab, 25, 0, 2000*GeV,
                    drlab, 25, 0, 5);

            histos2D[chan][label]["pt1_vs_pt2"] = bookHisto(chan + "_" + label + "_pt1_vs_pt2", label,
                    ptlab, 25, 0, 2000*GeV,
                    ptlab, 25, 0, 2000*GeV);

            // pt balance
            histos1D[chan][label]["pt1_minus_pt2"] = bookHisto(chan + "_" + label + "_pt1_minus_pt2", label,
                    ptlab, 25, -1000*GeV, 1000*GeV);
        }

        return;
    }


    void MC_BOOSTEDHBB::bookFourMomComp(const string& label1,
            const string& label2) {

        const string label = label1 + "_" + label2;

        string n;
        foreach (const string& chan, channels) {

            // dr and dpt
            histos1D[chan][label]["dr"] =
                bookHisto(chan + "_" + label + "_dr", drlab, label, 40, 0, 4.0);

            histos1D[chan][label]["dpt"] =
                bookHisto(chan + "_" + label + "_dpt", label,
                        ptlab, 25, -100*GeV, 100*GeV);


            // <pti> vs ptj
            n = "mean_" + label1 + "_pt_vs_" + label2 + "_pt";
            profiles1D[chan][label][n] =
                bookProfile(chan + "_" + n, label,
                        "$p_{T," + label2 + "} / GeV$", 25, 0, 500*GeV,
                        "$\\left< p_{T," + label1 + "} \\right>$");

            n = "mean_" + label2 + "_pt_vs_" + label1 + "_pt";
            profiles1D[chan][label][n] =
                bookProfile(chan + "_" + n, label,
                        "$p_{T," + label1 + "} / GeV$", 25, 0, 500*GeV,
                        "$\\left< p_{T," + label2 + "} \\right>$");


            // <dpt> vs pti
            n = "mean_" + label + "_dpt_vs_" + label1 + "_pt";
            profiles1D[chan][label][n] =
                bookProfile(chan + "_" + n, label,
                        "$p_{T," + label1 + "} / GeV$", 25, 0, 500*GeV,
                        "$\\left< p_{T," + label1 + "} - p_{T," + label2 + "} \\right>$");

            n = "mean_" + label + "_dpt_vs_" + label2 + "_pt";
            profiles1D[chan][label][n] =
                bookProfile(chan + "_" + n, label,
                        "$p_{T," + label2 + "} / GeV$", 25, 0, 500*GeV,
                        "$\\left< p_{T," + label1 + "} - p_{T," + label2 + "} \\right>$");



            // <dr> vs pti
            n = "mean_" + label + "_dr_vs_" + label1 + "_pt";
            profiles1D[chan][label][n] =
                bookProfile(chan + "_" + n, label,
                        "$p_{T," + label1 + "} / GeV$", 25, 0, 500*GeV,
                        "$\\left< \\Delta R(" + label1 + "," + label2 + "} \\right>$");

            n = "mean_" + label + "_dr_vs_" + label2 + "_pt";
            profiles1D[chan][label][n] =
                bookProfile(chan + "_" + n, label,
                        "$p_{T," + label2 + "} / GeV$", 25, 0, 500*GeV,
                        "$\\left< \\Delta R(" + label1 + "," + label2 + "} \\right>$");


            // <dr> vs dpt
            n = "mean_" + label + "_dr_vs_dpt";
            profiles1D[chan][label][n] =
                bookProfile(chan + "_" + n, label,
                        "$p_{T," + label1 + "}  - p_{T," + label2 + "} / GeV$", 25, -100*GeV, 100*GeV,
                        "$\\left< \\Delta R(" + label1 + "," + label2 + "} \\right>$");


            // <dpt> vs dr
            n = "mean_" + label + "_dpt_vs_dr";
            profiles1D[chan][label][n] =
                bookProfile(chan + "_" + n, label,
                        "$\\Delta R(" + label1 + "," + label2 + "}$", 40, 0, 4.0,
                        "\\left< $p_{T," + label1 + "}  - p_{T," + label2 + "} / GeV$ \\right>");



            // (dr, dpt)
            histos2D[chan][label]["dr_vs_dpt"] =
                bookHisto(chan + "_" + label + "_dr_vs_dpt", label,
                        ptlab, 25, -100*GeV, 100*GeV,
                        drlab, 40, 0, 4.0);

            // (pti, ptj)
            n = label1 + "_pt_vs_" + label2 + "_pt";
            histos2D[chan][label][n] = bookHisto(chan + "_" + n, label,
                    ptlab, 25, 0, 500*GeV,
                    ptlab, 25, 0, 500*GeV);
        }

        return;
    }



    void MC_BOOSTEDHBB::bookFourMomColl(const string& label) {
        bookFourMom(label);

        foreach (const string& chan, channels) {
            histos1D[chan][label]["n"] =
                bookHisto(chan + "_" + label + "_n", "multiplicity", "", 10, 0, 10);
        }

        return;
    }


    void MC_BOOSTEDHBB::bookBTagging(const string& label) {

        foreach (const string& chan, channels) {
            // book constituent fraction histograms
            histos1D[chan][label]["ConstitFracFromBHad"] =
                bookHisto(chan + "_" + label + "_ConstitFracFromBHad", label,
                        "jet constituent fraction from B hadron", 25, 0, 1.0);

            histos2D[chan][label]["ConstitFracFromBHad_vs_BHad_pt"] =
                bookHisto(chan + "_" + label + "_ConstitFracFromBHad", label,
                        ptlab, 25, 0, 2000*GeV,
                        "jet constituent fraction from B hadron", 25, 0, 1.0);

            histos2D[chan][label]["ConstitFracFromBHad_vs_Jet_pt"] =
                bookHisto(chan + "_" + label + "_ConstitFracFromBHad", label,
                        ptlab, 25, 0, 2000*GeV,
                        "jet constituent fraction from B hadron", 25, 0, 1.0);


            // book child fraction histograms
            histos1D[chan][label]["BHadChildFracInJet"] =
                bookHisto(chan + "_" + label + "_BHadChildFracInJet", label,
                        "B hadron child fraction in jet", 25, 0, 1.0);

            histos2D[chan][label]["BHadChildFracInJet_vs_BHad_pt"] =
                bookHisto(chan + "_" + label + "_BHadChildFracInJet_vs_BHad_pt", label,
                        ptlab, 25, 0, 2000*GeV,
                        "B hadron child fraction in jet", 25, 0, 1.0);

            histos2D[chan][label]["BHadChildFracInJet_vs_Jet_pt"] =
                bookHisto(chan + "_" + label + "_BHadChildFracInJet_vs_Jet_pt", label,
                        ptlab, 25, 0, 2000*GeV,
                        "B hadron child fraction in jet", 25, 0, 1.0);


            // book pt fraction histograms
            histos1D[chan][label]["PtFracFromBHad"] =
                bookHisto(chan + "_" + label + "_PtFracFromBHad", label,
                        "jet $p_T$ fraction from B hadron", 25, 0, 1.0);

            histos2D[chan][label]["PtFracFromBHad_vs_BHad_pt"] =
                bookHisto(chan + "_" + label + "_PtFracFromBHad", label,
                        ptlab, 25, 0, 2000*GeV,
                        "jet $p_T$ fraction from B hadron", 25, 0, 1.0);

            histos2D[chan][label]["PtFracFromBHad_vs_Jet_pt"] =
                bookHisto(chan + "_" + label + "_PtFracFromBHad", label,
                        ptlab, 25, 0, 2000*GeV,
                        "jet $p_T$ fraction from B hadron", 25, 0, 1.0);


            // book child pt fraction histograms
            histos1D[chan][label]["BHadChildPtFracInJet"] =
                bookHisto(chan + "_" + label + "_BHadChildPtFracInJet", label,
                        "B hadron child $p_T$ fraction in jet", 25, 0, 1.0);

            histos2D[chan][label]["BHadChildPtFracInJet_vs_BHad_pt"] =
                bookHisto(chan + "_" + label + "_BHadChildPtFracInJet_vs_BHad_pt", label,
                        ptlab, 25, 0, 2000*GeV,
                        "B hadron child $p_T$ fraction in jet", 25, 0, 1.0);

            histos2D[chan][label]["BHadChildPtFracInJet_vs_Jet_pt"] =
                bookHisto(chan + "_" + label + "_BHadChildPtFracInJet_vs_Jet_pt", label,
                        ptlab, 25, 0, 2000*GeV,
                        "B hadron child $p_T$ fraction in jet", 25, 0, 1.0);
        }


        return;
    }


    void MC_BOOSTEDHBB::fillFourMom(const string& channel, const string& label, const FourMomentum& p, double weight) {
        MSG_DEBUG("Filling " << label << " histograms");

        histos1D[channel][label]["pt"]->fill(p.pT(), weight);
        histos1D[channel][label]["eta"]->fill(p.eta(), weight);
        return;
    }

    void MC_BOOSTEDHBB::fillFourMomAllChannels(const string& label, const FourMomentum& p, double weight) {
        MSG_DEBUG("Filling " << label << " histograms");

        histos1DAllChannels[label]["pt"]->fill(p.pT(), weight);
        histos1DAllChannels[label]["eta"]->fill(p.eta(), weight);
        return;
    }


    void MC_BOOSTEDHBB::fillFourMomPair(const string& channel,
            const string& label1, const FourMomentum& p1,
            const string& label2, const FourMomentum& p2,
            double weight) {

        const string& label = label1 + "_" + label2;
        fillFourMom(channel, label, p1 + p2, weight);

        double dr = Rivet::deltaR(p1, p2);
        double pt = (p1 + p2).pT();

        histos1D[channel][label]["dr"]->fill(dr, weight);
        histos2D[channel][label]["dr_vs_pt"]->fill(pt, dr);
        histos2D[channel][label]["pt1_vs_pt2"]->fill(p2.pT(), p1.pT());
        histos1D[channel][label]["pt1_minus_pt2"]->fill(p1.pT() - p2.pT());

        return;
    }


    void MC_BOOSTEDHBB::fillFourMomComp(const string& channel,
            const string& label1, const FourMomentum& p1,
            const string& label2, const FourMomentum& p2,
            double weight) {

        const string label = label1 + "_" + label2;

        double dr = Rivet::deltaR(p1, p2);
        double pt1 = p1.pT();
        double pt2 = p2.pT();
        double dpt = pt1 - pt2;

        histos1D[channel][label]["dr"]->fill(dr, weight);
        histos1D[channel][label]["dpt"]->fill(dpt, weight);

        string n = "mean_" + label1 + "_pt_vs_" + label2 + "_pt";
        profiles1D[channel][label][n]->fill(pt2, pt1, weight);

        n = "mean_" + label2 + "_pt_vs_" + label1 + "_pt";
        profiles1D[channel][label][n]->fill(pt1, pt2, weight);

        n = "mean_" + label + "_dpt_vs_" + label1 + "_pt";
        profiles1D[channel][label][n]->fill(pt1, dpt, weight);

        n = "mean_" + label + "_dpt_vs_" + label2 + "_pt";
        profiles1D[channel][label][n]->fill(pt2, dpt, weight);

        n = "mean_" + label + "_dr_vs_" + label1 + "_pt";
        profiles1D[channel][label][n]->fill(pt1, dr, weight);

        n = "mean_" + label + "_dr_vs_" + label2 + "_pt";
        profiles1D[channel][label][n]->fill(pt2, dr, weight);

        n = "mean_" + label + "_dr_vs_dpt";
        profiles1D[channel][label][n]->fill(dpt, dr, weight);

        n = "mean_" + label + "_dpt_vs_dr";
        profiles1D[channel][label][n]->fill(dr, dpt, weight);

        histos2D[channel][label]["dr_vs_dpt"]->fill(dpt, dr, weight);

        n = label1 + "_pt_vs_" + label2 + "_pt";
        histos2D[channel][label][n]->fill(pt2, pt1);

        return;
    }


    template <class T>
    void MC_BOOSTEDHBB::fillFourMomColl(const string& channel,
            const string& label, const vector<T>& ps,
            double weight) {

        MSG_DEBUG("Filling " << ps.size() << " members of collection " << label);
        histos1D[channel][label]["n"]->fill(ps.size(), weight);

        foreach (const T& p, ps)
            fillFourMom(channel, label, p.mom(), weight);

        return;
    }


    void MC_BOOSTEDHBB::fillBTagging(const string& channel,
            const string& label, const Jet& jet, const Particle& bhad,
            double weight) {

        // TODO
        // visible and charged particle fractions?

        const Particles& commonParts = constsFromPart(jet, bhad);

        double constitFrac = double(commonParts.size())/jet.constituents().size();
        double childFrac = double(commonParts.size())/bhad.stableDescendants().size();

        double commonPartPt = 0.0;
        foreach (const Particle& p, commonParts)
            commonPartPt += p.pt();

        double constitPtFrac = commonPartPt/jet.pt();
        double childPtFrac = commonPartPt/bhad.pt();

        // fill constituent fraction histograms
        histos1D[channel][label]["ConstitFracFromBHad"]->fill(constitFrac, weight);

        histos2D[channel][label]["ConstitFracFromBHad_vs_BHad_pt"]->fill(bhad.pt(), constitFrac, weight);

        histos2D[channel][label]["ConstitFracFromBHad_vs_Jet_pt"]->fill(jet.pt(), constitFrac, weight);


        // fill child fraction histograms
        histos1D[channel][label]["BHadChildFracInJet"]->fill(childFrac, weight);

        histos2D[channel][label]["BHadChildFracInJet_vs_BHad_pt"]->fill(bhad.pt(), childFrac, weight);

        histos2D[channel][label]["BHadChildFracInJet_vs_Jet_pt"]->fill(jet.pt(), childFrac, weight);


        // fill pt fraction histograms
        histos1D[channel][label]["PtFracFromBHad"]->fill(constitPtFrac, weight);

        histos2D[channel][label]["PtFracFromBHad_vs_BHad_pt"]->fill(bhad.pt(), constitPtFrac, weight);

        histos2D[channel][label]["PtFracFromBHad_vs_Jet_pt"]->fill(jet.pt(), constitPtFrac, weight);


        // fill child pt fraction histograms
        histos1D[channel][label]["BHadChildPtFracInJet"]->fill(childPtFrac, weight);

        histos2D[channel][label]["BHadChildPtFracInJet_vs_BHad_pt"]->fill(bhad.pt(), childPtFrac, weight);

        histos2D[channel][label]["BHadChildPtFracInJet_vs_Jet_pt"]->fill(jet.pt(), childPtFrac, weight);


        return;
    }



    Jets MC_BOOSTEDHBB::bTagged(const Jets& js) {
        Jets bjs;
        foreach(const Jet& j, js) {
            if (j.bTagged())
                bjs.push_back(j);
        }

        return bjs;
    }


    void MC_BOOSTEDHBB::fillDRBHadHists(const Event& event,
            const Particle& bhad, const vector<JetCollection>& jetColls) {

        int jidx;
        double weight = event.weight();

        string name;
        double ptMin;
        foreach (const JetCollection& jetColl, jetColls) {
            name = jetColl.first;
            ptMin = jetColl.second;

            const Jets& jets =
                applyProjection<FastJets>(event, name).jetsByPt(ptMin);

            jidx = drMinJetIdx(bhad, jets);

            if (jidx >= 0) {
                fillFourMomComp(name, "DRBHad", bhad,
                        "Jet", jets[jidx], weight);
                fillFourMom(name, "DRBHad", bhad, weight);
                fillFourMom(name, "DRBHadJet", jets[jidx], weight);
            }

        }

        return;
    }


    void MC_BOOSTEDHBB::fillGABHadHists(const Event& event,
            const vector<JetCollection>& jetColls) {

        double weight = event.weight();

        string name;
        double ptMin;
        foreach (const JetCollection& jetColl, jetColls) {
            name = jetColl.first;
            ptMin = jetColl.second;

            const Jets& jets =
                applyProjection<FastJets>(event, name).jetsByPt(ptMin);

            cout << "about to loop over all " << name << " jets." << endl;
            foreach (const Jet& jet, jets) {
                cout << "found a jet with pt eta phi:" << endl
                    << jet.pt() << " " << jet.eta() << " " << jet.phi() << endl;

                //Note that jet.bTags() will return the b hadrons associated with that jet via ghost association.
                foreach (const Particle& bhad, jet.bTags()) {
                    fillFourMomComp(name, "GABHad", bhad,
                            "Jet", jet, weight);
                    fillFourMom(name, "GABHad", bhad, weight);
                    fillFourMom(name, "GABHadJet", jet, weight);
                }
            }

        }

        return;
    }

    void MC_BOOSTEDHBB::fillBHadAssociated(const Event& event,
            const vector<JetCollection>& jetColls) {

        double weight = event.weight();

        string name;
        double ptMin;
        foreach (const JetCollection& jetColl, jetColls) {
            name = jetColl.first;
            ptMin = jetColl.second;

            const Jets& jets =
                applyProjection<FastJets>(event, name).jetsByPt(ptMin);

            foreach (const Jet& jet, jets) {
                //Note that jet.bTags() will return the b hadrons associated with that jet via ghost association.
                foreach (const Particle& bhad, jet.bTags()) {
                    fillFourMom(name, "Associated-BHad", bhad, weight);
                }
            }

        }

        return;
    }


    const Particles MC_BOOSTEDHBB::constsFromPart(const Jet& jet, const Particle& bhad) const {
        const Particles& BChildren = bhad.stableDescendants();

        // TODO
        // I'm sure this is sloooooow.
        Particles parts;
        foreach (const Particle& p, jet.constituents()) {
            foreach (const Particle& child, BChildren) {
                if (p == child) {
                    parts.push_back(p);
                    break;
                }
            }
        }

        return parts;
    }




    //@}

} // Rivet
