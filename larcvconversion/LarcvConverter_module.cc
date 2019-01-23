////////////////////////////////////////////////////////////////////////
// Class:       LarcvConverter
// Plugin Type: analyzer (art v2_11_03)
// File:        larcv_converter_module.cc
//
// Generated at Tue Jan 22 11:26:24 2019 by Adams, Corey James using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//  Larcv includes:
#include "larcv/core/DataFormat/IOManager.h"

#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventVoxel2D.h"
#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"
#include "larcv/core/DataFormat/Voxel3DMeta.h"


// Larsoft includes:
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// larcv includes:


class LarcvConverter;


class LarcvConverter : public art::EDAnalyzer {
public:
  explicit LarcvConverter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LarcvConverter(LarcvConverter const &) = delete;
  LarcvConverter(LarcvConverter &&) = delete;
  LarcvConverter & operator = (LarcvConverter const &) = delete;
  LarcvConverter & operator = (LarcvConverter &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.
    void neutrino_slice(art::Event const & e, larcv::IOManager* io);

    void cluster_slice(art::Event const & e, larcv::IOManager *io);

    // void wire_slice(art::Event const & e, larcv::IOManager *io);

    void build_particle_map(art::Event const & e, larcv::IOManager* io);

    void initialize_meta();

    void reset_larcv_variables();

    std::vector< std::vector< int> > _particle_to_trackID;
    std::map< int, int > _trackID_to_particle;


    std::vector<larcv::ImageMeta> _image_meta_2d;
    larcv::Voxel3DMeta  _voxel_meta;

    // Private pointer to an io manager:
    larcv::IOManager * _io;


    // Variables from fcl:
    std::string larcv_out_name;
};


LarcvConverter::LarcvConverter(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  ,
  larcv_out_name(p.get<std::string>("larcv_out_name"))
 // More initializers here.
{}

void LarcvConverter::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  reset_larcv_variables();

  // Get the event ID information for this event:
  int run = e.eventAuxiliary().run();
  int subrun = e.eventAuxiliary().subRun();
  int event = e.eventAuxiliary().event();


  neutrino_slice(e, _io);



  // Save the event
  _io->set_id(run, subrun, event);
  _io->save_entry();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// This code needs to be in an intialize call somehow:
void LarcvConverter::beginJob(){

  _io = new IOManager(kWRITE);
  _io -> set_out_file("larcv_converted.root");

  if (!_io) {
    std::cout << "Must set io manager before initializing!" << std::endl;
    throw std::exception();
  }

  // init all modules:
  for (size_t n = 0; n < _modules.size(); n++) {
    _modules[n]->initialize();
  }

  _io->initialize();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// This code needs to be in an finalize call somehow:
void LarcvConverter::endJob(){

  _io->finalize();
  delete _io;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void LarcvConverter::reset_larcv_variables(){
  _particle_to_trackID.clear();
  _trackID_to_particle.clear();
}

void LarcvConverter::initialize_meta(){
  _voxel_meta = larcv::Voxel3DMeta();
  
  float voxel_3d_size = 1.0;
  _voxel_meta.set(-1000, -1000, 0, 1000, 1000, 2000, 2000/voxel_3d_size, 2000/voxel_3d_size, 2000/voxel_3d_size);

  _image_meta_2d.clear();

  // parameters for ImageMeta are (xmin, ymin, xmax, ymax, nx, ny, units)
  // We'll encode tick in y and wire in x.  Units will be centimeters
  // y (drift direction) goes from -200 to 200 for n_ticks * 2 + spacing
  // x (wire direction) goes from 0
  _max_tick = 4*n_ticks;

  // int _readout_length = 4492;
  // int _n_channels = 30720;

  _image_meta_2d.push_back(larcv::ImageMeta(
      0, 0, 30720, _max_tick, _max_tick / compression, 30720, 0, larcv::kUnitCM));
  _image_meta_2d.push_back(larcv::ImageMeta(
      0, 0, 30720, _max_tick, _max_tick / compression, 30720, 1, larcv::kUnitCM));
  _image_meta_2d.push_back(larcv::ImageMeta(
      0, 0, 30720, _max_tick, _max_tick / compression, 30720, 2, larcv::kUnitCM));

  return;
}






void LarcvConverter::build_particle_map(art::Event const & e, larcv::IOManager* io) {
  // This function makes the mapping between geant objects and larcv particles

  // It builds the list of particles in larcv, and populates the maps
  // _particle_to_trackID
  // _trackID_to_particle

  _particle_to_trackID.clear();
  _trackID_to_particle.clear();

  // Get the MCtracks and MCShowers

  std::string producer = "mcreco";
  art::InputTag tag(producer);
  auto const& mctracks = e.getValidHandle<std::vector<sim::MCTrack> >(tag);
  auto const& mcshowers = e.getValidHandle<std::vector<sim::MCShower> >(tag);

  // Get the EventParticle from larcv:
  auto event_part = (larcv::EventParticle*)io->get_data("particle", "duneseg");

  // std::cout << "Number of mctracks : " << mctracks->size() << std::endl;
  // std::cout << "Number of mcshowers: " << mcshowers->size() << std::endl;

  unsigned int id = 0;

  for (auto& track : *mctracks) {
    larcv::Particle part;

    part.id(id);
    part.mcst_index(          track.TrackID());

    part.track_id(            track.TrackID());
    part.pdg_code(            track.PdgCode());
    part.nu_interaction_type( track.Origin());
    part.creation_process(    track.Process());

    part.parent_track_id(     track.MotherTrackID());
    part.parent_pdg_code(     track.MotherPdgCode());
    part.ancestor_track_id(   track.AncestorTrackID());
    part.ancestor_pdg_code(   track.AncestorPdgCode());

    part.first_step( track.Start().Position().X(),
                     track.Start().Position().Y(),
                     track.Start().Position().Z(),
                     track.Start().Position().T());

    part.energy_init(track.Start().Momentum().E());
    part.momentum(   track.Start().Momentum().X(),
                     track.Start().Momentum().Y(),
                     track.Start().Momentum().Z());

    _particle_to_trackID.push_back(std::vector<int>());
    _particle_to_trackID.back().push_back(track.TrackID());
    _trackID_to_particle[track.TrackID()] = id;

    event_part->emplace_back(std::move(part));
    id++;
  }

  for (auto& shower : *mcshowers) {
    larcv::Particle part;

    part.id(id);
    part.mcst_index(          shower.TrackID());

    part.track_id(            shower.TrackID());
    part.pdg_code(            shower.PdgCode());
    part.nu_interaction_type( shower.Origin());
    part.creation_process(    shower.Process());

    part.parent_track_id(     shower.MotherTrackID());
    part.parent_pdg_code(     shower.MotherPdgCode());
    part.ancestor_track_id(   shower.AncestorTrackID());
    part.ancestor_pdg_code(   shower.AncestorPdgCode());

    part.first_step( shower.Start().Position().X(),
                     shower.Start().Position().Y(),
                     shower.Start().Position().Z(),
                     shower.Start().Position().T());

    part.energy_init(shower.Start().Momentum().E());
    part.momentum(   shower.Start().Momentum().X(),
                     shower.Start().Momentum().Y(),
                     shower.Start().Momentum().Z());

    _particle_to_trackID.push_back(std::vector<int>());
    for (auto& daughter_id : shower.DaughterTrackID()) {
      _particle_to_trackID.back().push_back(daughter_id);
      _trackID_to_particle[daughter_id] = id;
    }

    event_part->emplace_back(std::move(part));
    id++;
  }

  return;
}




void LarcvConverter::neutrino_slice(art::Event const & e, larcv::IOManager* io){

  std::string neutrino_producer = "generator";
  art::InputTag neutrino_tag(neutrino_producer);

  art::Handle<std::vector<simb::MCTruth> > mctruth;

  e->getByLabel(neutrino_tag, mctruth);

  auto truth = mctruth->at(0);
  auto neutrino = mctruth->at(0).GetNeutrino().Nu();

  auto event_particle  =
      (larcv::EventParticle*) io->get_data("particle", "duneneutrino");

  // Start by extracting the particle information:
  larcv::Particle neut_info;
  neut_info.id(0);

  // Info from MCNeutrino:
  neut_info.nu_interaction_type(truth.GetNeutrino().InteractionType());
  neut_info.nu_current_type(truth.GetNeutrino().CCNC());
  neut_info.track_id(neutrino.TrackId());
  neut_info.pdg_code(neutrino.PdgCode());
  neut_info.creation_process(neutrino.Process());
  neut_info.position(neutrino.Vx(),
                     neutrino.Vy(),
                     neutrino.Vz(),
                     neutrino.T());
  neut_info.momentum(neutrino.Px(),
                     neutrino.Py(),
                     neutrino.Pz());
  neut_info.energy_init(neutrino.E());

  event_particle->emplace_back(std::move(neut_info));

  for (size_t i = 0; i < truth.NParticles(); i ++){
    larcv::Particle particle;
    auto & larsoft_particle = truth.GetParticle(i);
    if (larsoft_particle.StatusCode() != 1){
      continue;
    }

    particle.id(i+1);
    particle.track_id(larsoft_particle.TrackId());
    particle.pdg_code(larsoft_particle.PdgCode());
    particle.parent_track_id(larsoft_particle.Mother());
    particle.creation_process(larsoft_particle.Process());

    particle.position(larsoft_particle.Vx(),
                       larsoft_particle.Vy(),
                       larsoft_particle.Vz(),
                       larsoft_particle.T());

    particle.momentum(larsoft_particle.Px(),
                       larsoft_particle.Py(),
                       larsoft_particle.Pz());

    particle.energy_init(larsoft_particle.E());
    event_particle->emplace_back(std::move(particle));
  }


  return;
}

void LarcvConverter::cluster_slice(art::Event const & e, larcv::IOManager *io){

  build_particle_map(e, io);


  // Get the simch data:
  std::string _simch_producer = "largeant";
  art::InputTag digit_tag(_simch_producer);
  auto const& simch =
      e.getValidHandle<std::vector<sim::SimChannel> >(digit_tag);


  // get the cluster3d objects:
  auto event_cluster3d =
      (larcv::EventClusterVoxel3D*)io->get_data("cluster3d", "duneseg");

  // This is to hold the clusters in 2d:
  auto event_cluster2d =
      (larcv::EventClusterPixel2D*)io->get_data("cluster2d", "duneseg");

  auto event_image2d = 
      (larcv::EventImage2D*)io->get_data("image2d", "dunewire");

  // This sets up placeholder images:
  std::vector<larcv::Image2D> images;
  for (size_t i = 0; i < _image_meta_2d.size(); i++)
    images.push_back(larcv::Image2D(_image_meta_2d.at(i)));


  // Now, loop over the sim channels, and add the depositions to the
  // correct voxels

  int n_particles = _particle_to_trackID.size();

  std::vector<larcv::ClusterPixel2D> _clusters_by_projection;
  _clusters_by_projection.resize(3);

  int i = 0;
  for (auto& cluster2dSet : _clusters_by_projection) {
    cluster2dSet.resize(n_particles + 1);
    cluster2dSet.meta(_image_meta_2d.at(i));
    i++;
  }

  // larcv::ClusterVoxel3D clusters3d;
  event_cluster3d->resize(n_particles + 1);
  event_cluster3d->meta(_voxel_meta);

  float _min_x(9999), _max_x(-9999);
  float _min_y(9999), _max_y(-9999);
  float _min_z(9999), _max_z(-9999);

  for (auto& ch : *simch) {
    int this_column = column(ch.Channel());
    int this_projection_id = projection_id(ch.Channel());

    for (auto& TDCIDE : ch.TDCIDEMap()) {
      auto& tdc = TDCIDE.first;
      auto& ides = TDCIDE.second;


      for (auto& ide : ides) {

        // if (tdc < 0 || tdc > n_ticks){
          // continue;
        // }
        if (ide.trackID == -1){
          continue;
        }

        // Add this ide to the proper particle in 3D:
        int this_particle = ide.trackID;
        int larcv_particle_id;

        if (_trackID_to_particle.find(this_particle) !=
            _trackID_to_particle.end()) {
          larcv_particle_id = _trackID_to_particle[this_particle];
        } else{
          larcv_particle_id = n_particles;
        }
        // if (tdc > 0 && tdc <= 3000){
          event_cluster3d->writeable_voxel_set(larcv_particle_id)
              .add(larcv::Voxel(_voxel_meta.id(ide.x, ide.y, ide.z), ide.energy));
        // }

        if (ide.x > _max_x) _max_x = ide.x;
        if (ide.y > _max_y) _max_y = ide.y;
        if (ide.z > _max_z) _max_z = ide.z;

        if (ide.x < _min_x) _min_x = ide.x;
        if (ide.y < _min_y) _min_y = ide.y;
        if (ide.z < _min_z) _min_z = ide.z;


        int tick = tdc;

        // if (tdc < 3000 && tdc > 0){
        // Write the energy to the voxel set:
        _clusters_by_projection.at(this_projection_id)
            .writeable_voxel_set(larcv_particle_id)
            .add(larcv::Voxel(
                _image_meta_2d.at(this_projection_id).index(tick / compression, this_column),
                ide.energy));
        // }

        // Write the energy to the wire image:
        auto index = images.at(this_projection_id).meta().index(tick / compression, this_column);
        float value = images.at(this_projection_id).pixel(tick / compression, this_column);
        images.at(this_projection_id).set_pixel(index, value + ide.energy);

      }
      //
    }
  }

  std::cout << "IDE Ranges: \n"
            << "\tx: (" <<_min_x << ", " << _max_x << ")\n"
            << "\ty: (" <<_min_y << ", " << _max_y << ")\n"
            << "\tz: (" <<_min_z << ", " << _max_z << ")\n";

  // for (auto& cluster_pix_2d : _clusters_by_projection) {
  //   event_cluster2d->emplace(std::move(cluster_pix_2d));
    // }
  // event_cluster3d->set(clusters3d, voxel_meta);
  //   std::cout << ch.TDCIDEMap().size() << std::endl;

  //   //  TDCIDE is std::pair<unsigned short, std::vector<sim::IDE> >
  //   //  TDCIDEMap is std::vector<TDCIDE>
  //   //  So, basically, it's a vector of pairs where the first element of
  //   the
  //   //  pair
  //   //  is the TDC, and the second element of the pair is the
  //   //  list of sim::IDEs that deposited there

  //   // // Loop over the digit and compress it:
  //   // for(size_t row = 0; row < 750; row ++){
  //   //   float val = digit.ADC(row*4);
  //   //   val += digit.ADC(row*4 + 1);
  //   //   val += digit.ADC(row*4 + 2);
  //   //   val += digit.ADC(row*4 + 3);
  //   //   // std::cout << "Setting at (" << column << ", " << row << ")" <<
  //   //   std::endl;
  //   //   _images.at(projection_id).set_pixel(column, row, val * 0.25);
  //   // }
  // }

  // // // Emplace the images:
  // // event_image2d -> emplace(std::move(_images));
  return

}










DEFINE_ART_MODULE(LarcvConverter)
