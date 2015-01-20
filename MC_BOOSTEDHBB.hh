// -*- C++ -*-
#ifndef RIVET_MC_BOOSTEDHBB_HH
#define RIVET_MC_BOOSTEDHBB_HH

#include "Rivet/Analysis.hh"

namespace Rivet {

    typedef pair<string, double> JetCollection;

    class MC_BOOSTEDHBB : public Analysis {
        public:
            /// Constructor
            MC_BOOSTEDHBB()
                : Analysis("MC_BOOSTEDHBB") { }

            /// Book histograms and initialise projections before the run
            void init();

            /// Perform the per-event analysis
            void analyze(const Event& event);

            /// Normalise histograms etc., after the run
            void finalize();

        private:

            vector<string> channels;
            vector<JetCollection> collections;
            map<string, map<string, map<string, Histo1DPtr> > > histos1D;
            map<string, map<string, Histo1DPtr> >  histos1DAllChannels;
            map<string, map<string, map<string, Histo2DPtr> > > histos2D;
            map<string, map<string, map<string, Profile1DPtr> > > profiles1D;


            void bookJetcollection(const string& collName);
            void bookChannel(const string& channel);

            Histo1DPtr bookHisto(const string& name, const string& title,
                    const string& xlabel, int nxbins, double xmin, double xmax);

            Histo2DPtr bookHisto(const string& name, const string& title,
                    const string& xlabel, int nxbins, double xmin, double xmax,
                    const string& ylabel, int nybins, double ymin, double ymax);

            Profile1DPtr bookProfile(const string& name, const string& title,
                    const string& xlabel, int nxbins, double xmin, double xmax,
                    const string& ylabel);


            void bookFourMom(const string& name);
            void bookFourMomAllChannels(const string& label);
            void bookFourMomPair(const string& name1, const string& name2);
            void bookFourMomComp(const string& name, const string& name2);
            void bookFourMomColl(const string& name);
            void bookBTagging(const string& name);

            void fillFourMomAllChannels(const string& label, const FourMomentum& p, double weight);

            void fillFourMom(const string& channel,
                    const string& name,
                    const FourMomentum& p,
                    double weight);

            void fillFourMomPair(const string& channel,
                    const string& name1,
                    const FourMomentum& p1,
                    const string& name2,
                    const FourMomentum& p2,
                    double weight);

            void fillFourMomComp(const string& channel,
                    const string& name1,
                    const FourMomentum& p1,
                    const string& name2,
                    const FourMomentum& p2,
                    double weight);

            template <class T>
            void fillFourMomColl(const string& channel,
                    const string& name,
                    const vector<T>& ps,
                    double weight);


            void fillBTagging(const string& channel,
                    const string& name,
                    const Jet& jet,
                    const Particle& bhad,
                    double weight);


            Jets bTagged(const Jets& js);

            void fillDRBHadHists(const Event& event,
                    const Particle& bhad,
                    const vector<JetCollection>& jetColls);

            void fillGABHadHists(const Event& event,
                    const vector<JetCollection>& jetColls);

            void fillBHadAssociated(const Event& event,
                    const vector<JetCollection>& jetColls);

            const Particles constsFromPart(const Jet& jet, const Particle& bhad) const;

    };


    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_BOOSTEDHBB);
}

#endif
