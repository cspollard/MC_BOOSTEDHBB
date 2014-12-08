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

        bookChannel("Rho05Min00Max10");
        bookChannel("Rho10Min00Max10");
        bookChannel("Rho20Min00Max10");
        bookChannel("Rho30Min00Max10");

        bookChannel("Rho05Min005Max10");
        bookChannel("Rho10Min005Max10");
        bookChannel("Rho20Min005Max10");
        bookChannel("Rho30Min005Max10");

        bookChannel("Rho05Min01Max10");
        bookChannel("Rho10Min01Max10");
        bookChannel("Rho20Min01Max10");
        bookChannel("Rho30Min01Max10");

        bookChannel("AKTTrack02");
        bookChannel("AKTTrack03");
        bookChannel("AKTTrack04");
        bookChannel("AKTCalo04");


        // calo jets constituents
        FinalState caloParts(-2.5, 2.5, 0.5*GeV);

        // track jets constituents
        ChargedFinalState trackParts(-2.5, 2.5, 0.5*GeV);


        // variable-R jet projections
        addProjection(FastJets(trackParts, aktVRPlugin(5*GeV, 0, 1)),
                "Rho05Min00Max10");
        addProjection(FastJets(trackParts, aktVRPlugin(10*GeV, 0, 1)),
                "Rho10Min00Max10");
        addProjection(FastJets(trackParts, aktVRPlugin(20*GeV, 0, 1)),
                "Rho20Min00Max10");
        addProjection(FastJets(trackParts, aktVRPlugin(30*GeV, 0, 1)),
                "Rho30Min00Max10");

        addProjection(FastJets(trackParts, aktVRPlugin(5*GeV, 0.05, 1)),
                "Rho05Min005Max10");
        addProjection(FastJets(trackParts, aktVRPlugin(10*GeV, 0.05, 1)),
                "Rho10Min005Max10");
        addProjection(FastJets(trackParts, aktVRPlugin(20*GeV, 0.05, 1)),
                "Rho20Min005Max10");
        addProjection(FastJets(trackParts, aktVRPlugin(30*GeV, 0.05, 1)),
                "Rho30Min005Max10");

        addProjection(FastJets(trackParts, aktVRPlugin(5*GeV, 0.1, 1)),
                "Rho05Min01Max10");
        addProjection(FastJets(trackParts, aktVRPlugin(10*GeV, 0.1, 1)),
                "Rho10Min01Max10");
        addProjection(FastJets(trackParts, aktVRPlugin(20*GeV, 0.1, 1)),
                "Rho20Min01Max10");
        addProjection(FastJets(trackParts, aktVRPlugin(30*GeV, 0.1, 1)),
                "Rho30Min01Max10");


        // conventional jet projections
        addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.2), "AKTTrack02");
        addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.3), "AKTTrack03");
        addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.4), "AKTTrack04");

        addProjection(FastJets(caloParts, FastJets::ANTIKT, 0.4), "AKTCalo04");


        // projection to find b-hadrons
        addProjection(HeavyHadrons(-2.5, 2.5, 5*GeV), "HeavyHadrons");

        //Now book histograms to plot
        bookFourMomComp("DRBHad", "Jet");
        bookFourMomComp("GABHad", "Jet");
        bookFourMom("DRBHad");
        bookFourMom("DRBHadJet");
        bookFourMom("GABHad");
        bookFourMom("GABHadJet");

        return;
    }


    /// Perform the per-event analysis
    void MC_BOOSTEDHBB::analyze(const Event& event) {
        const Particles& bhads =
            applyProjection<HeavyHadrons>(event, "HeavyHadrons").bHadrons();

        if (!bhads.size())
            vetoEvent;


        fillGABHadHists(event, channels);

        foreach (const Particle& bhad, bhads)
            fillDRBHadHists(event, bhad, channels);


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
            histos1D[chan][label]["m"] = bookHisto(chan + "_" + label + "_m", label, mlab, 25, 0, 1000*GeV);

            histos2D[chan][label]["m_vs_pt"] = bookHisto(chan + "_" + label + "_m_vs_pt", label,
                    ptlab, 25, 0, 2000*GeV,
                    mlab, 25, 0, 1000*GeV);
        }

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


    void MC_BOOSTEDHBB::fillFourMom(const string& channel, const string& label, const FourMomentum& p, double weight) {
        MSG_DEBUG("Filling " << label << " histograms");

        histos1D[channel][label]["pt"]->fill(p.pT(), weight);
        histos1D[channel][label]["eta"]->fill(p.eta(), weight);
        histos1D[channel][label]["m"]->fill(p.mass(), weight);
        histos2D[channel][label]["m_vs_pt"]->fill(p.pT(), p.mass(), weight);

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


    Jets MC_BOOSTEDHBB::bTagged(const Jets& js) {
        Jets bjs;
        foreach(const Jet& j, js) {
            if (j.bTagged())
                bjs.push_back(j);
        }

        return bjs;
    }


    void MC_BOOSTEDHBB::fillDRBHadHists(const Event& event,
            const Particle& bhad, const vector<string>& jetColls) {

        int jidx;
        double weight = event.weight();

        foreach (const string& jetColl, jetColls) {

            const Jets& jets =
                applyProjection<FastJets>(event, jetColl).jetsByPt(25*GeV);

            jidx = drMinJetIdx(bhad, jets);

            if (jidx >= 0) {
                fillFourMomComp(jetColl, "DRBHad", bhad,
                        "Jet", jets[jidx], weight);
                fillFourMom(jetColl, "DRBHad", bhad, weight);
                fillFourMom(jetColl, "DRBHadJet", jets[jidx], weight);
            }

        }

        return;
    }


    void MC_BOOSTEDHBB::fillGABHadHists(const Event& event,
            const vector<string>& jetColls) {

        double weight = event.weight();

        foreach (const string& jetColl, jetColls) {

            const Jets& jets =
                applyProjection<FastJets>(event, jetColl).jetsByPt(25*GeV);

            foreach (const Jet& jet, jets) {
                foreach (const Particle& bhad, jet.bTags()) {
                    fillFourMomComp(jetColl, "GABHad", bhad,
                            "Jet", jet, weight);
                    fillFourMom(jetColl, "GABHad", bhad, weight);
                    fillFourMom(jetColl, "GABHadJet", jet, weight);
                }
            }

        }

        return;
    }



    //@}

} // Rivet
