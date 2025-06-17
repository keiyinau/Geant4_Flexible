#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TRandom3.h>
#include <ROOT/RDataFrame.hxx>
auto Data_Filter(std::string TreeName, std::string filename, std::string Particle) {
    ROOT::RDataFrame df(TreeName, filename.c_str());
    auto filtered_df = df.Filter([Particle](const ROOT::VecOps::RVec<Char_t>& name) {
        return std::string(name.data()) == Particle;
    }, {"ParticleName"});
    return filtered_df;
}
void Data_Converter_Material() {
    ROOT::EnableImplicitMT();
    // Create a TFile to store data
    const std::string filename = "gamma_CsI_5000_127.root";
    const std::string Particle_Consideration[]={"e+","gamma"};


    for(const auto& Particle_ : Particle_Consideration) {
        TFile* file=new TFile(filename.c_str(), "READ");
        //Filter for Reference data
        auto ref_data=Data_Filter("Reference;1",filename,Particle_);
        std::vector<std::string> branches = ref_data.GetColumnNames();
        // Define a helper column for ranking stepID within (eventID, trackID)
        auto df_with_ranks = ref_data.Define("row_idx", "rdfentry_").Define("step_idx", [](const ROOT::VecOps::RVec<int>& s, const ROOT::VecOps::RVec<int>& e, const ROOT::VecOps::RVec<int>& t, ULong64_t i) {
            static std::map<std::pair<int, int>, std::pair<size_t, size_t>> m;
            if (m.empty() || i == 0) m.clear();
            auto k = std::make_pair(e[0], t[0]); auto& p = m[k];
            if (!p.first || s[0] < s[p.first - 1]) p.first = i + 1;
            if (!p.second || s[0] > s[p.second - 1]) p.second = i + 1;
            return p;
        }, {"StepID", "EventID", "TrackID", "row_idx"});
        df_with_ranks.Filter("row_idx + 1 == step_idx.first").Snapshot("InitialData", "Reference_"+Particle_+"_Initial.root", branches);
        df_with_ranks.Filter("row_idx + 1 == step_idx.second").Snapshot("FinalData", "Reference_"+Particle_+"_Final.root", branches);
        
        auto tracker_data=Data_Filter("Tracker;1",filename,Particle_);
        
        tracker_data.Snapshot("Tracker", "Tracker_"+Particle_+".root");
        file->Close();
    }

}
