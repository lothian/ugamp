#include "hamiltonian.h"
#include "globals.h"
#include <libiwl/iwl.h>
#include <libmints/wavefunction.h>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libtrans/integraltransform.h>
#include <libdpd/dpd.h>

#define ID(x) ints.DPD_ID(x)

namespace psi { namespace ugamp {

Hamiltonian::Hamiltonian(boost::shared_ptr<PSIO> psio, boost::shared_ptr<Wavefunction> ref,
std::vector<boost::shared_ptr<MOSpace> > spaces, bool ovov_only, bool full_virtual_space)
{
  Dimension frzvpi = ref->frzvpi();
  if(full_virtual_space)
    for(int h=0; h < ref->nirrep(); h++) frzvpi[h] = 0;

  nmo_ = ref->nmo();
  nfzc_ = ref->nfrzc();
  nfzv_ = 0;
  no_ = nv_ = 0;
  for(int i=0; i < ref->nirrep(); i++) {
    nfzv_ += frzvpi[i];
    no_ += ref->doccpi()[i] - ref->frzcpi()[i];
    nv_ += ref->nmopi()[i] - ref->doccpi()[i] - ref->frzvpi()[i];
  }
  nact_ = nmo_ - nfzc_ - nfzv_;

  int nact = nact_;
  int no = no_;
  int nv = nv_;
  ovov_only_ = ovov_only;

  outfile->Printf("\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t*       HAMILTONIAN      *\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\n");

  outfile->Printf("\tNMO    = %d\n", nmo_);
  outfile->Printf("\tNFZC   = %d\n", nfzc_);
  outfile->Printf("\tNFZV   = %d\n", nfzv_);
  outfile->Printf("\tNACT   = %d\n", nact_);
  outfile->Printf("\tNO     = %d\n", no_);
  outfile->Printf("\tNV     = %d\n", nv_);

  // Prepare Fock matrix in MO basis
  SharedMatrix Fa = ref->Fa()->clone();
  SharedMatrix Ca = ref->Ca();
  Fa->transform(Ca);

  int *map = init_int_array(nmo_); // Translates from Pitzer (including frozen docc) to QT
  reorder_qt((int *) ref->doccpi(), (int *) ref->soccpi(), (int *) ref->frzcpi(), (int *) frzvpi, 
             map, (int *) ref->nmopi(), ref->nirrep());

  fock_ = block_matrix(nact, nact);
  int mo_offset=0;
  for(int h=0; h < ref->nirrep(); h++) {
    int nmo = ref->nmopi()[h]; int nfv = frzvpi[h]; int nfc = ref->frzcpi()[h];
    for(int p=nfc; p < nmo-nfv; p++) {
      for(int q=nfc; q < nmo-nfv; q++) {
      int P = map[p+mo_offset]; int Q = map[q+mo_offset];
      fock_[P-nfzc_][Q-nfzc_] = Fa->get(h,p,q);
      }
    }
    mo_offset += nmo;
  }
  free(map);

  // Use reorder_qt() to generate a new mapping array w/o frozen core or virtual orbitals
  int *doccpi = init_int_array(ref->nirrep());
  int *nmopi = init_int_array(ref->nirrep());
  int *null = init_int_array(ref->nirrep());
  for(int h=0; h < ref->nirrep(); h++) {
    doccpi[h] = ref->doccpi()[h] - ref->frzcpi()[h];
    nmopi[h] = ref->nmopi()[h] - ref->frzcpi()[h] - frzvpi[h];
  }
  int *map2 = init_int_array(nact); // Translates from Pitzer (w/o frozen MOs) to QT
  reorder_qt(doccpi, (int *) ref->soccpi(), null, null, map2, nmopi, ref->nirrep());
  free(null); free(nmopi); free(doccpi);

  if(ovov_only) {
    IntegralTransform ints(ref, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly, IntegralTransform::QTOrder, 
                           full_virtual_space ? IntegralTransform::OccOnly : IntegralTransform::OccAndVir);
    ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    ints_ = init_4d_array(no, no, nv, nv);
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    dpdbuf4 K;
    dpd_set_default(ints.get_dpd_id());
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    for(int h=0; h < ref->nirrep(); h++) {
      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);
      for(int pq=0; pq < K.params->rowtot[h]; pq++) {
        int p = map2[ K.params->roworb[h][pq][0] ];
        int q = map2[ K.params->roworb[h][pq][1] ];
        for(int rs=0; rs < K.params->coltot[h]; rs++) {
          int r = map2[ K.params->colorb[h][rs][0] ];
          int s = map2[ K.params->colorb[h][rs][1] ];
          ints_[p][r][q][s] = K.matrix[h][pq][rs];
        }
      }
      global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, 1);
  }
  else {
    IntegralTransform ints(ref, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly, IntegralTransform::PitzerOrder, 
                           full_virtual_space ? IntegralTransform::OccOnly : IntegralTransform::OccAndVir);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    ints_ = init_4d_array(nact, nact, nact, nact);
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    dpdbuf4 K;
    dpd_set_default(ints.get_dpd_id());
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    for(int h=0; h < ref->nirrep(); h++) {
      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);
      for(int pq=0; pq < K.params->rowtot[h]; pq++) {
        int p = map2[ K.params->roworb[h][pq][0] ];
        int q = map2[ K.params->roworb[h][pq][1] ];
        for(int rs=0; rs < K.params->coltot[h]; rs++) {
          int r = map2[ K.params->colorb[h][rs][0] ];
          int s = map2[ K.params->colorb[h][rs][1] ];
          ints_[p][r][q][s] = K.matrix[h][pq][rs];
        }
      }
      global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, 1);
  }

  // L(pqrs) = 2<pq|rs> - <pq|sr>  
  if(ovov_only) {
    L_ = init_4d_array(no, no, nv, nv);
    for(int p=0; p < no; p++)
      for(int q=0; q < no; q++)
        for(int r=0; r < nv; r++)
          for(int s=0; s < nv; s++)
            L_[p][q][r][s] = 2*ints_[p][q][r][s] - ints_[p][q][s][r];
  }
  else {
    L_ = init_4d_array(nact, nact, nact, nact);
    for(int p=0; p < nact; p++)
      for(int q=0; q < nact; q++)
        for(int r=0; r < nact; r++)
          for(int s=0; s < nact; s++)
            L_[p][q][r][s] = 2*ints_[p][q][r][s] - ints_[p][q][s][r];
  }

}

Hamiltonian::~Hamiltonian()
{
  if(ovov_only_) {
    free_4d_array(ints_, no_, no_, nv_);
    free_4d_array(L_, no_, no_, nv_);
  }
  else {
    free_4d_array(ints_, nact_, nact_, nact_);
    free_4d_array(L_, nact_, nact_, nact_);
  }
  free_block(fock_); 
}

}} // namespace psi::ugamp
