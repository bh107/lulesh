//
//  LULESH Bohrium
//
//
//  Created by Steffen Karlsson, 12/05/15.
//  Copyright 2015 University of Copenhagen. All rights reserved.
//

#include <bxx/bohrium.hpp>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

using namespace bxx;

enum { VolumeError = -1, QStopError = -2 };

/* Progress */
#define LULESH_SHOW_PROGRESS 0

/* Boundary conditions */
#define XI_M        0x003
#define XI_M_SYMM   0x001
#define XI_M_FREE   0x002

#define XI_P        0x00c
#define XI_P_SYMM   0x004
#define XI_P_FREE   0x008

#define ETA_M       0x030
#define ETA_M_SYMM  0x010
#define ETA_M_FREE  0x020

#define ETA_P       0x0c0
#define ETA_P_SYMM  0x040
#define ETA_P_FREE  0x080

#define ZETA_M      0x300
#define ZETA_M_SYMM 0x100
#define ZETA_M_FREE 0x200

#define ZETA_P      0xc00
#define ZETA_P_SYMM 0x400
#define ZETA_P_FREE 0x800


/****************************************************/
/* Allow flexibility for arithmetic representations */
/****************************************************/
typedef float real4;
typedef double real8;
typedef long double real10;
typedef uint64_t Index_t;
typedef real8 Real_t;
typedef int Int_t;

inline real4  CBRT(real4  arg) { return cbrtf(arg) ; }
inline real8  CBRT(real8  arg) { return cbrt(arg) ; }
inline real10 CBRT(real10 arg) { return cbrtl(arg) ; }

inline real4  FABS(real4  arg) { return fabsf(arg); }
inline real8  FABS(real8  arg) { return fabs(arg); }
inline real10 FABS(real10 arg) { return fabsl(arg); }

inline real4  SQRT(real4  arg) { return sqrtf(arg) ; }
inline real8  SQRT(real8  arg) { return sqrt(arg) ; }
inline real10 SQRT(real10 arg) { return sqrtl(arg) ; }

/************************************************************/
/* Allow for flexible data layout experiments by separating */
/* array interface from underlying implementation.          */
/************************************************************/
struct Domain {

public:

    void AllocateNodalPersistent(size_t size) {
        m_x = zeros<Real_t>(size);
        m_y = zeros<Real_t>(size);
        m_z = zeros<Real_t>(size);

        m_xd = zeros<Real_t>(size);
        m_yd = zeros<Real_t>(size);
        m_zd = zeros<Real_t>(size);

        m_xdd = zeros<Real_t>(size);
        m_ydd = zeros<Real_t>(size);
        m_zdd = zeros<Real_t>(size);

        m_fx = zeros<Real_t>(size);
        m_fy = zeros<Real_t>(size);
        m_fz = zeros<Real_t>(size);

        m_nodalMass = zeros<Real_t>(size);
    }

    void AllocateElemPersistent(size_t size) {
        m_nodelist = zeros<Index_t>(8*size);

        m_lxim = zeros<Index_t>(size);
        m_lxip = zeros<Index_t>(size);
        m_letam = zeros<Index_t>(size);
        m_letap = zeros<Index_t>(size);
        m_lzetam = zeros<Index_t>(size);
        m_lzetap = zeros<Index_t>(size);

        m_elemBC = zeros<Int_t>(size);

        m_e = zeros<Real_t>(size);

        m_p = ones<Real_t>(size);
        m_q = ones<Real_t>(size);
        m_ql = zeros<Real_t>(size);
        m_qq = zeros<Real_t>(size);

        m_v = ones<Real_t>(size);
        m_volo = zeros<Real_t>(size);
        m_delv = zeros<Real_t>(size);
        m_vdov = zeros<Real_t>(size);

        m_arealg = zeros<Real_t>(size);

        m_ss = zeros<Real_t>(size);

        m_elemMass = zeros<Real_t>(size);
    }

    void AllocateElemTemporary(size_t size) {
        m_dxx = zeros<Real_t>(size);
        m_dyy = zeros<Real_t>(size);
        m_dzz = zeros<Real_t>(size);

        m_delv_xi = zeros<Real_t>(size);
        m_delv_eta = zeros<Real_t>(size);
        m_delv_zeta = zeros<Real_t>(size);

        m_delx_xi = zeros<Real_t>(size);
        m_delx_eta = zeros<Real_t>(size);
        m_delx_zeta = zeros<Real_t>(size);

        m_vnew = zeros<Real_t>(size);
    }

    void AllocateNodesets(size_t size) {
        m_symmX = zeros<Index_t>(size);
        m_symmY = zeros<Index_t>(size);
        m_symmZ = zeros<Index_t>(size);
    }

    void AllocateNodeElemIndexes() {
        Index_t numElem = this->numElem();
        Index_t numNode = this->numNode();

        /* set up node-centered indexing of elements */
        m_nodeElemCount = zeros<Index_t>(numNode);
        for (Index_t i=0; i<numElem; ++i) {
            for (Index_t j=0; j < 8; ++j) {
                ++m_nodeElemCount[scalar<Index_t>(m_nodelist[Index_t(8)*i+j])];
            }
        }

        m_nodeElemStart = zeros<Index_t>(numNode);
        m_nodeElemStart[0] = 0;

        for (Index_t i=1; i < numNode; ++i)
            m_nodeElemStart[i] = m_nodeElemStart[i-1]+m_nodeElemCount[i-1];

        m_nodeElemCornerList = zeros<Index_t>(m_nodeElemStart[numNode-1]+
                                              m_nodeElemStart[numNode-1]);

        m_nodeElemCount = zeros<Index_t>(numNode);

        for (Index_t i=0; i < numElem; ++i) {
            for (Index_t j=0; j < 8; ++j) {
                Index_t m = scalar<Index_t>(m_nodelist[Index_t(8)*i+j]);
                Index_t k = i*8+j;
                Index_t offset = scalar<Index_t>(m_nodeElemStart[m])
                                 + scalar<Index_t>(m_nodeElemCount[m]);
                m_nodeElemCornerList[offset] = k;
                ++m_nodeElemCount[m];
            }
        }

        for (Index_t i=0; i < m_nodeElemCornerList.len(); ++i) {
            Index_t clv = scalar<Index_t>(m_nodeElemCornerList[i]);
            if ((clv < 0) || (clv > numElem*8)) {
                fprintf(stderr,
                        "AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!\n");
                exit(1);
            }
        }
    }

    /* coordinates */
    multi_array<Real_t> m_x;
    multi_array<Real_t> m_y;
    multi_array<Real_t> m_z;

    /* velocities */
    multi_array<Real_t> m_xd;
    multi_array<Real_t> m_yd;
    multi_array<Real_t> m_zd;

    /* accelerations */
    multi_array<Real_t> m_xdd;
    multi_array<Real_t> m_ydd;
    multi_array<Real_t> m_zdd;

    /* forces */
    multi_array<Real_t> m_fx;
    multi_array<Real_t> m_fy;
    multi_array<Real_t> m_fz;

    /* mass */
    multi_array<Real_t> m_nodalMass;

    /* symmetry plane nodesets */
    multi_array<Index_t> m_symmX;
    multi_array<Index_t> m_symmY;
    multi_array<Index_t> m_symmZ;

    multi_array<Index_t> m_nodeElemCount;
    multi_array<Index_t> m_nodeElemStart;
    multi_array<Index_t> m_nodeElemCornerList;

    /* element-centered */
    multi_array<Index_t> m_nodelist;

    /* element connectivity across each face */
    multi_array<Index_t> m_lxim;
    multi_array<Index_t> m_lxip;
    multi_array<Index_t> m_letam;
    multi_array<Index_t> m_letap;
    multi_array<Index_t> m_lzetam;
    multi_array<Index_t> m_lzetap;

    /* symmetry/free-surface flags for each elem face */
    multi_array<Int_t> m_elemBC;

    /* principal strains -- temporary */
    multi_array<Real_t> m_dxx;
    multi_array<Real_t> m_dyy;
    multi_array<Real_t> m_dzz;

    /* velocity gradient -- temporary */
    multi_array<Real_t> m_delv_xi;
    multi_array<Real_t> m_delv_eta;
    multi_array<Real_t> m_delv_zeta;

    /* coordinate gradient -- temporary */
    multi_array<Real_t> m_delx_xi;
    multi_array<Real_t> m_delx_eta;
    multi_array<Real_t> m_delx_zeta;

    /* energy */
    multi_array<Real_t> m_e;

    /* pressure */
    multi_array<Real_t> m_p;
    /* q */
    multi_array<Real_t> m_q;
    /* linear term for q */
    multi_array<Real_t> m_ql;
    /* quadratic term for q */
    multi_array<Real_t> m_qq;

    /* relative volume */
    multi_array<Real_t> m_v;
    /* reference volume */
    multi_array<Real_t> m_volo;
    /* new relative volume -- temporary */
    multi_array<Real_t> m_vnew;
    /* m_vnew - m_v */
    multi_array<Real_t> m_delv;
    /* volume derivative over volume */
    multi_array<Real_t> m_vdov;

    /* characteristic length of an element */
    multi_array<Real_t> m_arealg;

    /* "sound speed" */
    multi_array<Real_t> m_ss;

    /* mass */
    multi_array<Real_t> m_elemMass;


    /* Node-centered */
    Real_t elem_x(Index_t idx)    { return scalar<Real_t>(m_x[idx]); }
    Real_t elem_y(Index_t idx)    { return scalar<Real_t>(m_y[idx]); }
    Real_t elem_z(Index_t idx)    { return scalar<Real_t>(m_z[idx]); }

    Real_t elem_xd(Index_t idx)   { return scalar<Real_t>(m_xd[idx]); }
    Real_t elem_yd(Index_t idx)   { return scalar<Real_t>(m_yd[idx]); }
    Real_t elem_zd(Index_t idx)   { return scalar<Real_t>(m_zd[idx]); }

    Real_t elem_xdd(Index_t idx)  { return scalar<Real_t>(m_xdd[idx]); }
    Real_t elem_ydd(Index_t idx)  { return scalar<Real_t>(m_ydd[idx]); }
    Real_t elem_zdd(Index_t idx)  { return scalar<Real_t>(m_zdd[idx]); }

    Real_t elem_fx(Index_t idx)   { return scalar<Real_t>(m_fx[idx]); }
    Real_t elem_fy(Index_t idx)   { return scalar<Real_t>(m_fy[idx]); }
    Real_t elem_fz(Index_t idx)   { return scalar<Real_t>(m_fz[idx]); }

    Real_t elem_nodalMass(Index_t idx) { return scalar<Real_t>(m_nodalMass[idx]); }

    Index_t elem_symmX(Index_t idx) { return scalar<Index_t>(m_symmX[idx]); }
    Index_t elem_symmY(Index_t idx) { return scalar<Index_t>(m_symmY[idx]); }
    Index_t elem_symmZ(Index_t idx) { return scalar<Index_t>(m_symmZ[idx]); }

    Index_t elem_nodeElemCount(Index_t idx)    { return scalar<Index_t>(m_nodeElemCount[idx]); }
    Index_t elem_nodeElemStart(Index_t idx)    { return scalar<Index_t>(m_nodeElemStart[idx]); }
    Index_t elem_nodeElemCornerList(Index_t i) { return scalar<Index_t>(m_nodeElemCornerList[i]); }

    /* Element-centered */
    // Index_t* nodelist(Index_t idx)         { return &m_nodelist[Index_t(8)*idx]; }

    Index_t  elem_lxim(Index_t idx)   { return scalar<Index_t>(m_lxim[idx]); }
    Index_t  elem_lxip(Index_t idx)   { return scalar<Index_t>(m_lxip[idx]); }
    Index_t  elem_letam(Index_t idx)  { return scalar<Index_t>(m_letam[idx]); }
    Index_t  elem_letap(Index_t idx)  { return scalar<Index_t>(m_letap[idx]); }
    Index_t  elem_lzetam(Index_t idx) { return scalar<Index_t>(m_lzetam[idx]); }
    Index_t  elem_lzetap(Index_t idx) { return scalar<Index_t>(m_lzetap[idx]); }

    Int_t  elem_elemBC(Index_t idx) { return scalar<Int_t>(m_elemBC[idx]); }

    Real_t elem_dxx(Index_t idx)  { return scalar<Real_t>(m_dxx[idx]); }
    Real_t elem_dyy(Index_t idx)  { return scalar<Real_t>(m_dyy[idx]); }
    Real_t elem_dzz(Index_t idx)  { return scalar<Real_t>(m_dzz[idx]); }

    Real_t elem_delv_xi(Index_t idx)    { return scalar<Real_t>(m_delv_xi[idx]); }
    Real_t elem_delv_eta(Index_t idx)   { return scalar<Real_t>(m_delv_eta[idx]); }
    Real_t elem_delv_zeta(Index_t idx)  { return scalar<Real_t>(m_delv_zeta[idx]); }

    Real_t elem_delx_xi(Index_t idx)    { return scalar<Real_t>(m_delx_xi[idx]); }
    Real_t elem_delx_eta(Index_t idx)   { return scalar<Real_t>(m_delx_eta[idx]); }
    Real_t elem_delx_zeta(Index_t idx)  { return scalar<Real_t>(m_delx_zeta[idx]); }

    Real_t elem_e(Index_t idx)          { return scalar<Real_t>(m_e[idx]); }

    Real_t elem_p(Index_t idx)          { return scalar<Real_t>(m_p[idx]); }
    Real_t elem_q(Index_t idx)          { return scalar<Real_t>(m_q[idx]); }
    Real_t elem_ql(Index_t idx)         { return scalar<Real_t>(m_ql[idx]); }
    Real_t elem_qq(Index_t idx)         { return scalar<Real_t>(m_qq[idx]); }

    Real_t elem_v(Index_t idx)          { return scalar<Real_t>(m_v[idx]); }
    Real_t elem_volo(Index_t idx)       { return scalar<Real_t>(m_volo[idx]); }
    Real_t elem_vnew(Index_t idx)       { return scalar<Real_t>(m_vnew[idx]); }
    Real_t elem_delv(Index_t idx)       { return scalar<Real_t>(m_delv[idx]); }
    Real_t elem_vdov(Index_t idx)       { return scalar<Real_t>(m_vdov[idx]); }

    Real_t elem_arealg(Index_t idx) { return scalar<Real_t>(m_arealg[idx]); }

    Real_t elem_ss(Index_t idx) { return scalar<Real_t>(m_ss[idx]); }

    Real_t elem_elemMass(Index_t idx) { return scalar<Real_t>(m_elemMass[idx]); }

    /* Params */
    Real_t& dtfixed()              { return m_dtfixed; }
    Real_t& time()                 { return m_time; }
    Real_t& deltatime()            { return m_deltatime; }
    Real_t& deltatimemultlb()      { return m_deltatimemultlb; }
    Real_t& deltatimemultub()      { return m_deltatimemultub; }
    Real_t& stoptime()             { return m_stoptime; }

    Real_t& u_cut()                { return m_u_cut; }
    Real_t& hgcoef()               { return m_hgcoef; }
    Real_t& qstop()                { return m_qstop; }
    Real_t& monoq_max_slope()      { return m_monoq_max_slope; }
    Real_t& monoq_limiter_mult()   { return m_monoq_limiter_mult; }
    Real_t& e_cut()                { return m_e_cut; }
    Real_t& p_cut()                { return m_p_cut; }
    Real_t& ss4o3()                { return m_ss4o3; }
    Real_t& q_cut()                { return m_q_cut; }
    Real_t& v_cut()                { return m_v_cut; }
    Real_t& qlc_monoq()            { return m_qlc_monoq; }
    Real_t& qqc_monoq()            { return m_qqc_monoq; }
    Real_t& qqc()                  { return m_qqc; }
    Real_t& eosvmax()              { return m_eosvmax; }
    Real_t& eosvmin()              { return m_eosvmin; }
    Real_t& pmin()                 { return m_pmin; }
    Real_t& emin()                 { return m_emin; }
    Real_t& dvovmax()              { return m_dvovmax; }
    Real_t& refdens()              { return m_refdens; }

    Real_t& dtcourant()            { return m_dtcourant; }
    Real_t& dthydro()              { return m_dthydro; }
    Real_t& dtmax()                { return m_dtmax; }

    Int_t&  cycle()                { return m_cycle; }

    Index_t&  sizeX()              { return m_sizeX; }
    Index_t&  sizeY()              { return m_sizeY; }
    Index_t&  sizeZ()              { return m_sizeZ; }
    Index_t&  numElem()            { return m_numElem; }
    Index_t&  numNode()            { return m_numNode; }

private:

    /* fixed time increment */
    Real_t m_dtfixed;
    /* current time */
    Real_t m_time;
    /* variable time increment */
    Real_t m_deltatime;
    Real_t m_deltatimemultlb;
    Real_t m_deltatimemultub;
    /* end time for simulation */
    Real_t m_stoptime;

    /* velocity tolerance */
    Real_t m_u_cut;
    /* hourglass control */
    Real_t m_hgcoef;
    /* excessive q indicator */
    Real_t m_qstop;
    Real_t m_monoq_max_slope;
    Real_t m_monoq_limiter_mult;
    /* energy tolerance */
    Real_t m_e_cut;
    /* pressure tolerance */
    Real_t m_p_cut;
    Real_t m_ss4o3;
    /* q tolerance */
    Real_t m_q_cut;
    /* relative volume tolerance */
    Real_t m_v_cut;
    /* linear term coef for q */
    Real_t m_qlc_monoq;
    /* quadratic term coef for q */
    Real_t m_qqc_monoq;
    Real_t m_qqc;
    Real_t m_eosvmax;
    Real_t m_eosvmin;
    /* pressure floor */
    Real_t m_pmin;
    /* energy floor */
    Real_t m_emin;
    /* maximum allowable volume change */
    Real_t m_dvovmax;
    /* reference density */
    Real_t m_refdens;

    /* courant constraint */
    Real_t m_dtcourant;
    /* volume change constraint */
    Real_t m_dthydro;
    /* maximum allowable time increment */
    Real_t m_dtmax;

    /* iteration count for simulation */
    Int_t m_cycle;

    /* X,Y,Z extent of this block */
    Index_t m_sizeX;
    Index_t m_sizeY;
    Index_t m_sizeZ;

    /* Elements/Nodes in this domain */
    Index_t m_numElem;
    Index_t m_numNode;
} domain;

template <typename T>
T *Allocate(size_t size)
{
	return static_cast<T *>(malloc(sizeof(T)*size)) ;
}

template <typename T>
void Release(T **ptr)
{
	if (*ptr != NULL) {
		free(*ptr) ;
		*ptr = NULL ;
	}
}

static inline
Real_t CalcElemVolume(const Real_t x0, const Real_t x1,
                      const Real_t x2, const Real_t x3,
                      const Real_t x4, const Real_t x5,
                      const Real_t x6, const Real_t x7,
                      const Real_t y0, const Real_t y1,
                      const Real_t y2, const Real_t y3,
                      const Real_t y4, const Real_t y5,
                      const Real_t y6, const Real_t y7,
                      const Real_t z0, const Real_t z1,
                      const Real_t z2, const Real_t z3,
                      const Real_t z4, const Real_t z5,
                      const Real_t z6, const Real_t z7) {
    Real_t twelveth = Real_t(1.0)/Real_t(12.0);

    Real_t dx61 = x6 - x1;
    Real_t dy61 = y6 - y1;
    Real_t dz61 = z6 - z1;

    Real_t dx70 = x7 - x0;
    Real_t dy70 = y7 - y0;
    Real_t dz70 = z7 - z0;

    Real_t dx63 = x6 - x3;
    Real_t dy63 = y6 - y3;
    Real_t dz63 = z6 - z3;

    Real_t dx20 = x2 - x0;
    Real_t dy20 = y2 - y0;
    Real_t dz20 = z2 - z0;

    Real_t dx50 = x5 - x0;
    Real_t dy50 = y5 - y0;
    Real_t dz50 = z5 - z0;

    Real_t dx64 = x6 - x4;
    Real_t dy64 = y6 - y4;
    Real_t dz64 = z6 - z4;

    Real_t dx31 = x3 - x1;
    Real_t dy31 = y3 - y1;
    Real_t dz31 = z3 - z1;

    Real_t dx72 = x7 - x2;
    Real_t dy72 = y7 - y2;
    Real_t dz72 = z7 - z2;

    Real_t dx43 = x4 - x3;
    Real_t dy43 = y4 - y3;
    Real_t dz43 = z4 - z3;

    Real_t dx57 = x5 - x7;
    Real_t dy57 = y5 - y7;
    Real_t dz57 = z5 - z7;

    Real_t dx14 = x1 - x4;
    Real_t dy14 = y1 - y4;
    Real_t dz14 = z1 - z4;

    Real_t dx25 = x2 - x5;
    Real_t dy25 = y2 - y5;
    Real_t dz25 = z2 - z5;

#define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) \
((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

    Real_t volume =
    TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,
                   dy31 + dy72, dy63, dy20,
                   dz31 + dz72, dz63, dz20) +
    TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,
                   dy43 + dy57, dy64, dy70,
                   dz43 + dz57, dz64, dz70) +
    TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,
                   dy14 + dy25, dy61, dy50,
                   dz14 + dz25, dz61, dz50);

#undef TRIPLE_PRODUCT

    volume *= twelveth;

    return volume;
}

static inline
Real_t CalcElemVolume(const Real_t x[8], const Real_t y[8], const Real_t z[8])
{
    return CalcElemVolume(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
                          y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
                          z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
}

static inline
void CalcPositionForNodes(const Real_t dt) {
    domain.m_x += domain.m_xd * dt;
    domain.m_y += domain.m_yd * dt;
    domain.m_z += domain.m_zd * dt;
}

static inline
void CalcVelocityForNodes(const Real_t dt, const Real_t u_cut) {
    domain.m_xd(domain.m_xd + domain.m_xdd * dt);
    domain.m_xd(as<Real_t>(abs(domain.m_xd) >= u_cut)*domain.m_xd);

    domain.m_yd(domain.m_yd + domain.m_ydd * dt);
    domain.m_yd(as<Real_t>(abs(domain.m_yd) >= u_cut)*domain.m_yd);

    domain.m_zd(domain.m_zd + domain.m_zdd * dt);
    domain.m_zd(as<Real_t>(abs(domain.m_zd) >= u_cut)*domain.m_zd);
}

static inline
void ApplyAccelerationBoundaryConditionsForNodes() {
    //TODO: Scatter
    Real_t size = domain.sizeX();
    domain.m_xdd[_(0, -size, size+1)] = Real_t(0.0);

    Index_t sliceSize = domain.m_xdd.len() - (size*size) - (size+1);
    multi_array<Real_t> symmY;
    symmY = view_as(domain.m_ydd[_(0, sliceSize - 1)], sliceSize / size, size);
    symmY[_(0, sliceSize, (size+1)*(size+1))][_(0, -1)] = Real_t(0.0);

    domain.m_zdd[_(0, size*(size+2), 1)] = Real_t(0.0);

    // printf("m_xdd: %d, m_symmX: %d", domain.m_xdd.len(), domain.m_symmX.len());
    // std::cout << "Values: " <<  domain.m_symmX << std::endl;
    //
    // scatter(domain.m_xdd, zeros<Real_t>(domain.m_xdd.len()), domain.m_symmX);
    // scatter(domain.m_ydd, zeros<Real_t>(domain.m_ydd.len()), domain.m_symmY);
    // scatter(domain.m_zdd, zeros<Real_t>(domain.m_zdd.len()), domain.m_symmZ);
}

static inline
void CalcAccelerationForNodes() {
    printf("CalcAccelerationForNodes\n");
    domain.m_xdd(domain.m_fx / domain.m_nodalMass);
    domain.m_ydd(domain.m_fy / domain.m_nodalMass);
    domain.m_zdd(domain.m_fz / domain.m_nodalMass);

    printf("fx: ");
    for( Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode )
        printf("%f, ", scalar<Real_t>(domain.m_fx[gnode]));
    printf("\n");
    printf("fy: ");
    for( Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode )
        printf("%f, ", scalar<Real_t>(domain.m_fy[gnode]));
    printf("\n");
    printf("fz: ");
    for( Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode )
        printf("%f, ", scalar<Real_t>(domain.m_fz[gnode]));
    printf("\n");

    printf("m_xdd: ");
    for( Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode )
        printf("%f, ", scalar<Real_t>(domain.m_xdd[gnode]));
    printf("\n");
    printf("m_ydd: ");
    for( Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode )
        printf("%f, ", scalar<Real_t>(domain.m_ydd[gnode]));
    printf("\n");
    printf("m_zdd: ");
    for( Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode )
        printf("%f, ", scalar<Real_t>(domain.m_zdd[gnode]));
    printf("\n");

    printf("m_nodalMass: ");
    for( Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode )
        printf("%f, ", scalar<Real_t>(domain.m_nodalMass[gnode]));
    printf("\n");
}

static inline
void InitStressTermsForElems(multi_array<Real_t>& sigxx,
                             multi_array<Real_t>& sigyy,
                             multi_array<Real_t>& sigzz) {
    sigxx(sigyy(sigzz(((Real_t)-1) * (domain.m_p - domain.m_q))));
}

static inline
Real_t CalcElemShapeFunctionDerivatives(const Real_t* const x,
									    const Real_t* const y,
									    const Real_t* const z,
									    Real_t b[][8]) {
    const Real_t x0 = x[0] ;   const Real_t x1 = x[1] ;
	const Real_t x2 = x[2] ;   const Real_t x3 = x[3] ;
	const Real_t x4 = x[4] ;   const Real_t x5 = x[5] ;
	const Real_t x6 = x[6] ;   const Real_t x7 = x[7] ;

	const Real_t y0 = y[0] ;   const Real_t y1 = y[1] ;
	const Real_t y2 = y[2] ;   const Real_t y3 = y[3] ;
	const Real_t y4 = y[4] ;   const Real_t y5 = y[5] ;
	const Real_t y6 = y[6] ;   const Real_t y7 = y[7] ;

	const Real_t z0 = z[0] ;   const Real_t z1 = z[1] ;
	const Real_t z2 = z[2] ;   const Real_t z3 = z[3] ;
	const Real_t z4 = z[4] ;   const Real_t z5 = z[5] ;
	const Real_t z6 = z[6] ;   const Real_t z7 = z[7] ;

	Real_t fjxxi, fjxet, fjxze;
	Real_t fjyxi, fjyet, fjyze;
	Real_t fjzxi, fjzet, fjzze;
	Real_t cjxxi, cjxet, cjxze;
	Real_t cjyxi, cjyet, cjyze;
	Real_t cjzxi, cjzet, cjzze;

	fjxxi = .125 * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) );
	fjxet = .125 * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) );
	fjxze = .125 * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) );

	fjyxi = .125 * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) );
	fjyet = .125 * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) );
	fjyze = .125 * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) );

	fjzxi = .125 * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) );
	fjzet = .125 * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) );
	fjzze = .125 * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) );

	/* compute cofactors */
	cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
	cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
	cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);

	cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
	cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
	cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);

	cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
	cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
	cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);

	/* calculate partials :
     this need only be done for l = 0,1,2,3   since , by symmetry ,
     (6,7,4,5) = - (0,1,2,3) .
	 */
	b[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
	b[0][1] =      cjxxi  -  cjxet  -  cjxze;
	b[0][2] =      cjxxi  +  cjxet  -  cjxze;
	b[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
	b[0][4] = -b[0][2];
	b[0][5] = -b[0][3];
	b[0][6] = -b[0][0];
	b[0][7] = -b[0][1];

	b[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
	b[1][1] =      cjyxi  -  cjyet  -  cjyze;
	b[1][2] =      cjyxi  +  cjyet  -  cjyze;
	b[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
	b[1][4] = -b[1][2];
	b[1][5] = -b[1][3];
	b[1][6] = -b[1][0];
	b[1][7] = -b[1][1];

	b[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
	b[2][1] =      cjzxi  -  cjzet  -  cjzze;
	b[2][2] =      cjzxi  +  cjzet  -  cjzze;
	b[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
	b[2][4] = -b[2][2];
	b[2][5] = -b[2][3];
	b[2][6] = -b[2][0];
	b[2][7] = -b[2][1];

	/* calculate jacobian determinant (volume) */
	return Real_t(8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
}

static inline
void SumElemFaceNormal(Real_t *normalX0, Real_t *normalY0, Real_t *normalZ0,
                       Real_t *normalX1, Real_t *normalY1, Real_t *normalZ1,
                       Real_t *normalX2, Real_t *normalY2, Real_t *normalZ2,
                       Real_t *normalX3, Real_t *normalY3, Real_t *normalZ3,
                       const Real_t x0, const Real_t y0, const Real_t z0,
                       const Real_t x1, const Real_t y1, const Real_t z1,
                       const Real_t x2, const Real_t y2, const Real_t z2,
                       const Real_t x3, const Real_t y3, const Real_t z3) {
	Real_t bisectX0 = Real_t(0.5) * (x3 + x2 - x1 - x0);
	Real_t bisectY0 = Real_t(0.5) * (y3 + y2 - y1 - y0);
	Real_t bisectZ0 = Real_t(0.5) * (z3 + z2 - z1 - z0);
	Real_t bisectX1 = Real_t(0.5) * (x2 + x1 - x3 - x0);
	Real_t bisectY1 = Real_t(0.5) * (y2 + y1 - y3 - y0);
	Real_t bisectZ1 = Real_t(0.5) * (z2 + z1 - z3 - z0);
	Real_t areaX = Real_t(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
	Real_t areaY = Real_t(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
	Real_t areaZ = Real_t(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

	*normalX0 += areaX;
	*normalX1 += areaX;
	*normalX2 += areaX;
	*normalX3 += areaX;

	*normalY0 += areaY;
	*normalY1 += areaY;
	*normalY2 += areaY;
	*normalY3 += areaY;

	*normalZ0 += areaZ;
	*normalZ1 += areaZ;
	*normalZ2 += areaZ;
	*normalZ3 += areaZ;
}

static inline
void CalcElemNodeNormals(Real_t pfx[8],
                         Real_t pfy[8],
                         Real_t pfz[8],
                         const Real_t x[8],
                         const Real_t y[8],
                         const Real_t z[8]) {
	for (Index_t i = 0 ; i < 8 ; ++i) {
		pfx[i] = Real_t(0.0);
		pfy[i] = Real_t(0.0);
		pfz[i] = Real_t(0.0);
	}
	/* evaluate face one: nodes 0, 1, 2, 3 */
	SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
					  &pfx[1], &pfy[1], &pfz[1],
					  &pfx[2], &pfy[2], &pfz[2],
					  &pfx[3], &pfy[3], &pfz[3],
					  x[0], y[0], z[0], x[1], y[1], z[1],
					  x[2], y[2], z[2], x[3], y[3], z[3]);
	/* evaluate face two: nodes 0, 4, 5, 1 */
	SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
					  &pfx[4], &pfy[4], &pfz[4],
					  &pfx[5], &pfy[5], &pfz[5],
					  &pfx[1], &pfy[1], &pfz[1],
					  x[0], y[0], z[0], x[4], y[4], z[4],
					  x[5], y[5], z[5], x[1], y[1], z[1]);
	/* evaluate face three: nodes 1, 5, 6, 2 */
	SumElemFaceNormal(&pfx[1], &pfy[1], &pfz[1],
					  &pfx[5], &pfy[5], &pfz[5],
					  &pfx[6], &pfy[6], &pfz[6],
					  &pfx[2], &pfy[2], &pfz[2],
					  x[1], y[1], z[1], x[5], y[5], z[5],
					  x[6], y[6], z[6], x[2], y[2], z[2]);
	/* evaluate face four: nodes 2, 6, 7, 3 */
	SumElemFaceNormal(&pfx[2], &pfy[2], &pfz[2],
					  &pfx[6], &pfy[6], &pfz[6],
					  &pfx[7], &pfy[7], &pfz[7],
					  &pfx[3], &pfy[3], &pfz[3],
					  x[2], y[2], z[2], x[6], y[6], z[6],
					  x[7], y[7], z[7], x[3], y[3], z[3]);
	/* evaluate face five: nodes 3, 7, 4, 0 */
	SumElemFaceNormal(&pfx[3], &pfy[3], &pfz[3],
					  &pfx[7], &pfy[7], &pfz[7],
					  &pfx[4], &pfy[4], &pfz[4],
					  &pfx[0], &pfy[0], &pfz[0],
					  x[3], y[3], z[3], x[7], y[7], z[7],
					  x[4], y[4], z[4], x[0], y[0], z[0]);
	/* evaluate face six: nodes 4, 7, 6, 5 */
	SumElemFaceNormal(&pfx[4], &pfy[4], &pfz[4],
					  &pfx[7], &pfy[7], &pfz[7],
					  &pfx[6], &pfy[6], &pfz[6],
					  &pfx[5], &pfy[5], &pfz[5],
					  x[4], y[4], z[4], x[7], y[7], z[7],
					  x[6], y[6], z[6], x[5], y[5], z[5]);
}

static inline
void IntegrateStressForElems(multi_array<Real_t>& sigxx,
                             multi_array<Real_t>& sigyy,
                             multi_array<Real_t>& sigzz,
                             multi_array<Real_t>& determ,
                             Index_t numElem) {
    Index_t numElem8 = numElem * 8 ;
    Real_t *fx_elem = Allocate<Real_t>(numElem8) ;
    Real_t *fy_elem = Allocate<Real_t>(numElem8) ;
    Real_t *fz_elem = Allocate<Real_t>(numElem8) ;

//TODO: Optimize
    for(Index_t k=0 ; k<numElem ; ++k) {
        Real_t B[3][8];
		Real_t x_local[8];
		Real_t y_local[8];
		Real_t z_local[8];

        for( Index_t lnode=0 ; lnode<8 ; ++lnode ) {
            Index_t gnode = scalar<Index_t>(domain.m_nodelist[Index_t(8)*k+lnode]);
            x_local[lnode] = scalar<Real_t>(domain.m_x[gnode]);
			y_local[lnode] = scalar<Real_t>(domain.m_y[gnode]);
			z_local[lnode] = scalar<Real_t>(domain.m_z[gnode]);
		}
        determ[k] = CalcElemShapeFunctionDerivatives(x_local, y_local, z_local, B);

        CalcElemNodeNormals(B[0], B[1], B[2], x_local, y_local, z_local);

        Index_t k8 = k*8;
        Real_t stress_xx = scalar<Real_t>(sigxx[k]);
        Real_t stress_yy = scalar<Real_t>(sigyy[k]);
        Real_t stress_zz = scalar<Real_t>(sigzz[k]);

        fx_elem[k8+0] = -( stress_xx * B[0][0] );
        fx_elem[k8+1] = -( stress_xx * B[0][1] );
    	fx_elem[k8+2] = -( stress_xx * B[0][2] );
    	fx_elem[k8+3] = -( stress_xx * B[0][3] );
    	fx_elem[k8+4] = -( stress_xx * B[0][4] );
    	fx_elem[k8+5] = -( stress_xx * B[0][5] );
    	fx_elem[k8+6] = -( stress_xx * B[0][6] );
    	fx_elem[k8+7] = -( stress_xx * B[0][7] );

    	fy_elem[k8+0] = -( stress_yy * B[1][0] );
    	fy_elem[k8+1] = -( stress_yy * B[1][1] );
    	fy_elem[k8+2] = -( stress_yy * B[1][2] );
    	fy_elem[k8+3] = -( stress_yy * B[1][3] );
    	fy_elem[k8+4] = -( stress_yy * B[1][4] );
    	fy_elem[k8+5] = -( stress_yy * B[1][5] );
    	fy_elem[k8+6] = -( stress_yy * B[1][6] );
    	fy_elem[k8+7] = -( stress_yy * B[1][7] );

    	fz_elem[k8+0] = -( stress_zz * B[2][0] );
    	fz_elem[k8+1] = -( stress_zz * B[2][1] );
    	fz_elem[k8+2] = -( stress_zz * B[2][2] );
    	fz_elem[k8+3] = -( stress_zz * B[2][3] );
    	fz_elem[k8+4] = -( stress_zz * B[2][4] );
    	fz_elem[k8+5] = -( stress_zz * B[2][5] );
    	fz_elem[k8+6] = -( stress_zz * B[2][6] );
    	fz_elem[k8+7] = -( stress_zz * B[2][7] );
    }

    for(Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode) {
        Index_t count = scalar<Index_t>(domain.m_nodeElemCount[gnode]);
        Index_t start = scalar<Index_t>(domain.m_nodeElemStart[gnode]);
        Real_t fx = Real_t(0.0);
        Real_t fy = Real_t(0.0);
        Real_t fz = Real_t(0.0);
        for (Index_t i=0 ; i < count ; ++i) {
            Index_t elem = scalar<Index_t>(domain.m_nodeElemCornerList[start+i]);
            fx += fx_elem[elem];
            fy += fy_elem[elem];
            fz += fz_elem[elem];
        }
        domain.m_fx[gnode] = fx;
        domain.m_fy[gnode] = fy;
        domain.m_fz[gnode] = fz;
    }

    printf("m_fx:");
	for (Index_t i = 0 ; i < domain.numNode() ; ++i)
		printf("%f," , scalar<Real_t>(domain.m_fx[i]));
	printf("\nm_fy:");
	for (Index_t i = 0 ; i < domain.numNode() ; ++i)
		printf("%f," , scalar<Real_t>(domain.m_fy[i]));
    printf("\nm_fz:");
	for (Index_t i = 0 ; i < domain.numNode() ; ++i)
		printf("%f," , scalar<Real_t>(domain.m_fz[i]));
	printf("\n");

    Release(&fz_elem);
	Release(&fy_elem);
	Release(&fx_elem);
}

static inline
void CollectDomainNodesToElemNodes(Index_t i,
                                   Real_t elemX[8],
                                   Real_t elemY[8],
                                   Real_t elemZ[8]) {
    Index_t nd0i = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+0]);
    Index_t nd1i = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+1]);
    Index_t nd2i = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+2]);
    Index_t nd3i = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+3]);
    Index_t nd4i = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+4]);
    Index_t nd5i = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+5]);
    Index_t nd6i = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+6]);
    Index_t nd7i = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+7]);

    elemX[0] = scalar<Real_t>(domain.m_x[nd0i]);
    elemX[1] = scalar<Real_t>(domain.m_x[nd1i]);
    elemX[2] = scalar<Real_t>(domain.m_x[nd2i]);
    elemX[3] = scalar<Real_t>(domain.m_x[nd3i]);
    elemX[4] = scalar<Real_t>(domain.m_x[nd4i]);
    elemX[5] = scalar<Real_t>(domain.m_x[nd5i]);
    elemX[6] = scalar<Real_t>(domain.m_x[nd6i]);
    elemX[7] = scalar<Real_t>(domain.m_x[nd7i]);

    elemY[0] = scalar<Real_t>(domain.m_y[nd0i]);
    elemY[1] = scalar<Real_t>(domain.m_y[nd1i]);
    elemY[2] = scalar<Real_t>(domain.m_y[nd2i]);
    elemY[3] = scalar<Real_t>(domain.m_y[nd3i]);
    elemY[4] = scalar<Real_t>(domain.m_y[nd4i]);
    elemY[5] = scalar<Real_t>(domain.m_y[nd5i]);
    elemY[6] = scalar<Real_t>(domain.m_y[nd6i]);
    elemY[7] = scalar<Real_t>(domain.m_y[nd7i]);

    elemZ[0] = scalar<Real_t>(domain.m_z[nd0i]);
    elemZ[1] = scalar<Real_t>(domain.m_z[nd1i]);
    elemZ[2] = scalar<Real_t>(domain.m_z[nd2i]);
    elemZ[3] = scalar<Real_t>(domain.m_z[nd3i]);
    elemZ[4] = scalar<Real_t>(domain.m_z[nd4i]);
    elemZ[5] = scalar<Real_t>(domain.m_z[nd5i]);
    elemZ[6] = scalar<Real_t>(domain.m_z[nd6i]);
    elemZ[7] = scalar<Real_t>(domain.m_z[nd7i]);
}

static inline
void CalcElemFBHourglassForce(Real_t *xd, Real_t *yd, Real_t *zd,  Real_t *hourgam0,
                              Real_t *hourgam1, Real_t *hourgam2, Real_t *hourgam3,
                              Real_t *hourgam4, Real_t *hourgam5, Real_t *hourgam6,
                              Real_t *hourgam7, Real_t coefficient,
                              Real_t *hgfx, Real_t *hgfy, Real_t *hgfz ) {
	Index_t i00=0;
	Index_t i01=1;
	Index_t i02=2;
	Index_t i03=3;

	Real_t h00 =
	hourgam0[i00] * xd[0] + hourgam1[i00] * xd[1] +
	hourgam2[i00] * xd[2] + hourgam3[i00] * xd[3] +
	hourgam4[i00] * xd[4] + hourgam5[i00] * xd[5] +
	hourgam6[i00] * xd[6] + hourgam7[i00] * xd[7];

	Real_t h01 =
	hourgam0[i01] * xd[0] + hourgam1[i01] * xd[1] +
	hourgam2[i01] * xd[2] + hourgam3[i01] * xd[3] +
	hourgam4[i01] * xd[4] + hourgam5[i01] * xd[5] +
	hourgam6[i01] * xd[6] + hourgam7[i01] * xd[7];

	Real_t h02 =
	hourgam0[i02] * xd[0] + hourgam1[i02] * xd[1]+
	hourgam2[i02] * xd[2] + hourgam3[i02] * xd[3]+
	hourgam4[i02] * xd[4] + hourgam5[i02] * xd[5]+
	hourgam6[i02] * xd[6] + hourgam7[i02] * xd[7];

	Real_t h03 =
	hourgam0[i03] * xd[0] + hourgam1[i03] * xd[1] +
	hourgam2[i03] * xd[2] + hourgam3[i03] * xd[3] +
	hourgam4[i03] * xd[4] + hourgam5[i03] * xd[5] +
	hourgam6[i03] * xd[6] + hourgam7[i03] * xd[7];

	hgfx[0] = coefficient *
	(hourgam0[i00] * h00 + hourgam0[i01] * h01 +
	 hourgam0[i02] * h02 + hourgam0[i03] * h03);

	hgfx[1] = coefficient *
	(hourgam1[i00] * h00 + hourgam1[i01] * h01 +
	 hourgam1[i02] * h02 + hourgam1[i03] * h03);

	hgfx[2] = coefficient *
	(hourgam2[i00] * h00 + hourgam2[i01] * h01 +
	 hourgam2[i02] * h02 + hourgam2[i03] * h03);

	hgfx[3] = coefficient *
	(hourgam3[i00] * h00 + hourgam3[i01] * h01 +
	 hourgam3[i02] * h02 + hourgam3[i03] * h03);

	hgfx[4] = coefficient *
	(hourgam4[i00] * h00 + hourgam4[i01] * h01 +
	 hourgam4[i02] * h02 + hourgam4[i03] * h03);

	hgfx[5] = coefficient *
	(hourgam5[i00] * h00 + hourgam5[i01] * h01 +
	 hourgam5[i02] * h02 + hourgam5[i03] * h03);

	hgfx[6] = coefficient *
	(hourgam6[i00] * h00 + hourgam6[i01] * h01 +
	 hourgam6[i02] * h02 + hourgam6[i03] * h03);

	hgfx[7] = coefficient *
	(hourgam7[i00] * h00 + hourgam7[i01] * h01 +
	 hourgam7[i02] * h02 + hourgam7[i03] * h03);

	h00 =
	hourgam0[i00] * yd[0] + hourgam1[i00] * yd[1] +
	hourgam2[i00] * yd[2] + hourgam3[i00] * yd[3] +
	hourgam4[i00] * yd[4] + hourgam5[i00] * yd[5] +
	hourgam6[i00] * yd[6] + hourgam7[i00] * yd[7];

	h01 =
	hourgam0[i01] * yd[0] + hourgam1[i01] * yd[1] +
	hourgam2[i01] * yd[2] + hourgam3[i01] * yd[3] +
	hourgam4[i01] * yd[4] + hourgam5[i01] * yd[5] +
	hourgam6[i01] * yd[6] + hourgam7[i01] * yd[7];

	h02 =
	hourgam0[i02] * yd[0] + hourgam1[i02] * yd[1]+
	hourgam2[i02] * yd[2] + hourgam3[i02] * yd[3]+
	hourgam4[i02] * yd[4] + hourgam5[i02] * yd[5]+
	hourgam6[i02] * yd[6] + hourgam7[i02] * yd[7];

	h03 =
	hourgam0[i03] * yd[0] + hourgam1[i03] * yd[1] +
	hourgam2[i03] * yd[2] + hourgam3[i03] * yd[3] +
	hourgam4[i03] * yd[4] + hourgam5[i03] * yd[5] +
	hourgam6[i03] * yd[6] + hourgam7[i03] * yd[7];


	hgfy[0] = coefficient *
	(hourgam0[i00] * h00 + hourgam0[i01] * h01 +
	 hourgam0[i02] * h02 + hourgam0[i03] * h03);

	hgfy[1] = coefficient *
	(hourgam1[i00] * h00 + hourgam1[i01] * h01 +
	 hourgam1[i02] * h02 + hourgam1[i03] * h03);

	hgfy[2] = coefficient *
	(hourgam2[i00] * h00 + hourgam2[i01] * h01 +
	 hourgam2[i02] * h02 + hourgam2[i03] * h03);

	hgfy[3] = coefficient *
	(hourgam3[i00] * h00 + hourgam3[i01] * h01 +
	 hourgam3[i02] * h02 + hourgam3[i03] * h03);

	hgfy[4] = coefficient *
	(hourgam4[i00] * h00 + hourgam4[i01] * h01 +
	 hourgam4[i02] * h02 + hourgam4[i03] * h03);

	hgfy[5] = coefficient *
	(hourgam5[i00] * h00 + hourgam5[i01] * h01 +
	 hourgam5[i02] * h02 + hourgam5[i03] * h03);

	hgfy[6] = coefficient *
	(hourgam6[i00] * h00 + hourgam6[i01] * h01 +
	 hourgam6[i02] * h02 + hourgam6[i03] * h03);

	hgfy[7] = coefficient *
	(hourgam7[i00] * h00 + hourgam7[i01] * h01 +
	 hourgam7[i02] * h02 + hourgam7[i03] * h03);

	h00 =
	hourgam0[i00] * zd[0] + hourgam1[i00] * zd[1] +
	hourgam2[i00] * zd[2] + hourgam3[i00] * zd[3] +
	hourgam4[i00] * zd[4] + hourgam5[i00] * zd[5] +
	hourgam6[i00] * zd[6] + hourgam7[i00] * zd[7];

	h01 =
	hourgam0[i01] * zd[0] + hourgam1[i01] * zd[1] +
	hourgam2[i01] * zd[2] + hourgam3[i01] * zd[3] +
	hourgam4[i01] * zd[4] + hourgam5[i01] * zd[5] +
	hourgam6[i01] * zd[6] + hourgam7[i01] * zd[7];

	h02 =
	hourgam0[i02] * zd[0] + hourgam1[i02] * zd[1]+
	hourgam2[i02] * zd[2] + hourgam3[i02] * zd[3]+
	hourgam4[i02] * zd[4] + hourgam5[i02] * zd[5]+
	hourgam6[i02] * zd[6] + hourgam7[i02] * zd[7];

	h03 =
	hourgam0[i03] * zd[0] + hourgam1[i03] * zd[1] +
	hourgam2[i03] * zd[2] + hourgam3[i03] * zd[3] +
	hourgam4[i03] * zd[4] + hourgam5[i03] * zd[5] +
	hourgam6[i03] * zd[6] + hourgam7[i03] * zd[7];


	hgfz[0] = coefficient *
	(hourgam0[i00] * h00 + hourgam0[i01] * h01 +
	 hourgam0[i02] * h02 + hourgam0[i03] * h03);

	hgfz[1] = coefficient *
	(hourgam1[i00] * h00 + hourgam1[i01] * h01 +
	 hourgam1[i02] * h02 + hourgam1[i03] * h03);

	hgfz[2] = coefficient *
	(hourgam2[i00] * h00 + hourgam2[i01] * h01 +
	 hourgam2[i02] * h02 + hourgam2[i03] * h03);

	hgfz[3] = coefficient *
	(hourgam3[i00] * h00 + hourgam3[i01] * h01 +
	 hourgam3[i02] * h02 + hourgam3[i03] * h03);

	hgfz[4] = coefficient *
	(hourgam4[i00] * h00 + hourgam4[i01] * h01 +
	 hourgam4[i02] * h02 + hourgam4[i03] * h03);

	hgfz[5] = coefficient *
	(hourgam5[i00] * h00 + hourgam5[i01] * h01 +
	 hourgam5[i02] * h02 + hourgam5[i03] * h03);

	hgfz[6] = coefficient *
	(hourgam6[i00] * h00 + hourgam6[i01] * h01 +
	 hourgam6[i02] * h02 + hourgam6[i03] * h03);

	hgfz[7] = coefficient *
	(hourgam7[i00] * h00 + hourgam7[i01] * h01 +
	 hourgam7[i02] * h02 + hourgam7[i03] * h03);
}

static inline
void CalcFBHourglassForceForElems(multi_array<Real_t>& determ,
								  Real_t *x8n, Real_t *y8n, Real_t *z8n,
								  Real_t *dvdx, Real_t *dvdy, Real_t *dvdz,
								  Real_t hourg) {
//TODO: Optimize
    // printf("CalcFBHourglassForceForElems\n");
    Index_t numElem = domain.numElem();
    Index_t numElem8 = numElem * 8;
    Real_t *fx_elem = Allocate<Real_t>(numElem8) ;
	Real_t *fy_elem = Allocate<Real_t>(numElem8) ;
	Real_t *fz_elem = Allocate<Real_t>(numElem8) ;

    Real_t  gamma[4][8];

    gamma[0][0] = Real_t( 1.);
    gamma[0][1] = Real_t( 1.);
    gamma[0][2] = Real_t(-1.);
    gamma[0][3] = Real_t(-1.);
    gamma[0][4] = Real_t(-1.);
    gamma[0][5] = Real_t(-1.);
    gamma[0][6] = Real_t( 1.);
    gamma[0][7] = Real_t( 1.);
    gamma[1][0] = Real_t( 1.);
    gamma[1][1] = Real_t(-1.);
    gamma[1][2] = Real_t(-1.);
    gamma[1][3] = Real_t( 1.);
    gamma[1][4] = Real_t(-1.);
    gamma[1][5] = Real_t( 1.);
    gamma[1][6] = Real_t( 1.);
    gamma[1][7] = Real_t(-1.);
    gamma[2][0] = Real_t( 1.);
    gamma[2][1] = Real_t(-1.);
    gamma[2][2] = Real_t( 1.);
    gamma[2][3] = Real_t(-1.);
    gamma[2][4] = Real_t( 1.);
    gamma[2][5] = Real_t(-1.);
    gamma[2][6] = Real_t( 1.);
    gamma[2][7] = Real_t(-1.);
    gamma[3][0] = Real_t(-1.);
    gamma[3][1] = Real_t( 1.);
    gamma[3][2] = Real_t(-1.);
    gamma[3][3] = Real_t( 1.);
    gamma[3][4] = Real_t( 1.);
    gamma[3][5] = Real_t(-1.);
    gamma[3][6] = Real_t( 1.);
    gamma[3][7] = Real_t(-1.);

    /*************************************************/
    /*    compute the hourglass modes */


    for(Index_t i2=0; i2<numElem; ++i2) {
        Real_t *fx_local, *fy_local, *fz_local;
        Real_t hgfx[8], hgfy[8], hgfz[8] ;

        Real_t coefficient;

        Real_t hourgam0[4], hourgam1[4], hourgam2[4], hourgam3[4] ;
        Real_t hourgam4[4], hourgam5[4], hourgam6[4], hourgam7[4];
        Real_t xd1[8], yd1[8], zd1[8] ;

        Index_t i3=8*i2;
        Real_t volinv=Real_t(1.0)/scalar<Real_t>(determ[i2]);
        Real_t ss1, mass1, volume13 ;

        for(Index_t i1=0;i1<4;++i1) {
            Real_t hourmodx =
            x8n[i3] * gamma[i1][0] + x8n[i3+1] * gamma[i1][1] +
            x8n[i3+2] * gamma[i1][2] + x8n[i3+3] * gamma[i1][3] +
            x8n[i3+4] * gamma[i1][4] + x8n[i3+5] * gamma[i1][5] +
            x8n[i3+6] * gamma[i1][6] + x8n[i3+7] * gamma[i1][7];

            Real_t hourmody =
            y8n[i3] * gamma[i1][0] + y8n[i3+1] * gamma[i1][1] +
            y8n[i3+2] * gamma[i1][2] + y8n[i3+3] * gamma[i1][3] +
            y8n[i3+4] * gamma[i1][4] + y8n[i3+5] * gamma[i1][5] +
            y8n[i3+6] * gamma[i1][6] + y8n[i3+7] * gamma[i1][7];

            Real_t hourmodz =
            z8n[i3] * gamma[i1][0] + z8n[i3+1] * gamma[i1][1] +
            z8n[i3+2] * gamma[i1][2] + z8n[i3+3] * gamma[i1][3] +
            z8n[i3+4] * gamma[i1][4] + z8n[i3+5] * gamma[i1][5] +
            z8n[i3+6] * gamma[i1][6] + z8n[i3+7] * gamma[i1][7];

            hourgam0[i1] = gamma[i1][0] -  volinv*(dvdx[i3  ] * hourmodx +
                                                   dvdy[i3  ] * hourmody +
                                                   dvdz[i3  ] * hourmodz );

            hourgam1[i1] = gamma[i1][1] -  volinv*(dvdx[i3+1] * hourmodx +
                                                   dvdy[i3+1] * hourmody +
                                                   dvdz[i3+1] * hourmodz );

            hourgam2[i1] = gamma[i1][2] -  volinv*(dvdx[i3+2] * hourmodx +
                                                   dvdy[i3+2] * hourmody +
                                                   dvdz[i3+2] * hourmodz );

            hourgam3[i1] = gamma[i1][3] -  volinv*(dvdx[i3+3] * hourmodx +
                                                   dvdy[i3+3] * hourmody +
                                                   dvdz[i3+3] * hourmodz );

            hourgam4[i1] = gamma[i1][4] -  volinv*(dvdx[i3+4] * hourmodx +
                                                   dvdy[i3+4] * hourmody +
                                                   dvdz[i3+4] * hourmodz );

            hourgam5[i1] = gamma[i1][5] -  volinv*(dvdx[i3+5] * hourmodx +
                                                   dvdy[i3+5] * hourmody +
                                                   dvdz[i3+5] * hourmodz );

            hourgam6[i1] = gamma[i1][6] -  volinv*(dvdx[i3+6] * hourmodx +
                                                   dvdy[i3+6] * hourmody +
                                                   dvdz[i3+6] * hourmodz );

            hourgam7[i1] = gamma[i1][7] -  volinv*(dvdx[i3+7] * hourmodx +
                                                   dvdy[i3+7] * hourmody +
                                                   dvdz[i3+7] * hourmodz );

        }

        /* compute forces */
        /* store forces into h arrays (force arrays) */

        ss1 = scalar<Real_t>(domain.m_ss[i2]);
        mass1 = scalar<Real_t>(domain.m_elemMass[i2]);
        volume13 = CBRT(scalar<Real_t>(determ[i2]));

        Index_t n0si2 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i2+0]);
        Index_t n1si2 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i2+1]);
        Index_t n2si2 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i2+2]);
        Index_t n3si2 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i2+3]);
        Index_t n4si2 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i2+4]);
        Index_t n5si2 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i2+5]);
        Index_t n6si2 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i2+6]);
        Index_t n7si2 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i2+7]);

        xd1[0] = scalar<Real_t>(domain.m_xd[n0si2]);
        xd1[1] = scalar<Real_t>(domain.m_xd[n1si2]);
        xd1[2] = scalar<Real_t>(domain.m_xd[n2si2]);
        xd1[3] = scalar<Real_t>(domain.m_xd[n3si2]);
        xd1[4] = scalar<Real_t>(domain.m_xd[n4si2]);
        xd1[5] = scalar<Real_t>(domain.m_xd[n5si2]);
        xd1[6] = scalar<Real_t>(domain.m_xd[n6si2]);
        xd1[7] = scalar<Real_t>(domain.m_xd[n7si2]);

        yd1[0] = scalar<Real_t>(domain.m_yd[n0si2]);
        yd1[1] = scalar<Real_t>(domain.m_yd[n1si2]);
        yd1[2] = scalar<Real_t>(domain.m_yd[n2si2]);
        yd1[3] = scalar<Real_t>(domain.m_yd[n3si2]);
        yd1[4] = scalar<Real_t>(domain.m_yd[n4si2]);
        yd1[5] = scalar<Real_t>(domain.m_yd[n5si2]);
        yd1[6] = scalar<Real_t>(domain.m_yd[n6si2]);
        yd1[7] = scalar<Real_t>(domain.m_yd[n7si2]);

        zd1[0] = scalar<Real_t>(domain.m_zd[n0si2]);
        zd1[1] = scalar<Real_t>(domain.m_zd[n1si2]);
        zd1[2] = scalar<Real_t>(domain.m_zd[n2si2]);
        zd1[3] = scalar<Real_t>(domain.m_zd[n3si2]);
        zd1[4] = scalar<Real_t>(domain.m_zd[n4si2]);
        zd1[5] = scalar<Real_t>(domain.m_zd[n5si2]);
        zd1[6] = scalar<Real_t>(domain.m_zd[n6si2]);
        zd1[7] = scalar<Real_t>(domain.m_zd[n7si2]);

        coefficient = - hourg * Real_t(0.01) * ss1 * mass1 / volume13;

        CalcElemFBHourglassForce(xd1,yd1,zd1,
                                 hourgam0,hourgam1,hourgam2,hourgam3,
                                 hourgam4,hourgam5,hourgam6,hourgam7,
                                 coefficient, hgfx, hgfy, hgfz);

        fx_local = &fx_elem[i3];
        fx_local[0] = hgfx[0];
		fx_local[1] = hgfx[1];
		fx_local[2] = hgfx[2];
		fx_local[3] = hgfx[3];
		fx_local[4] = hgfx[4];
		fx_local[5] = hgfx[5];
		fx_local[6] = hgfx[6];
		fx_local[7] = hgfx[7];

        fy_local = &fy_elem[i3] ;
		fy_local[0] = hgfy[0];
		fy_local[1] = hgfy[1];
		fy_local[2] = hgfy[2];
		fy_local[3] = hgfy[3];
		fy_local[4] = hgfy[4];
		fy_local[5] = hgfy[5];
		fy_local[6] = hgfy[6];
		fy_local[7] = hgfy[7];

		fz_local = &fz_elem[i3] ;
		fz_local[0] = hgfz[0];
		fz_local[1] = hgfz[1];
		fz_local[2] = hgfz[2];
		fz_local[3] = hgfz[3];
		fz_local[4] = hgfz[4];
		fz_local[5] = hgfz[5];
		fz_local[6] = hgfz[6];
		fz_local[7] = hgfz[7];
    }

    for( Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode ) {
        Index_t count = scalar<Index_t>(domain.m_nodeElemCount[gnode]);
        Index_t start = scalar<Index_t>(domain.m_nodeElemStart[gnode]);
        Real_t fx = Real_t(0.0);
        Real_t fy = Real_t(0.0);
        Real_t fz = Real_t(0.0);
        for (Index_t i=0 ; i < count ; ++i) {
            Index_t elem = scalar<Index_t>(domain.m_nodeElemCornerList[start+i]);
            fx += fx_elem[elem];
            fy += fy_elem[elem];
            fz += fz_elem[elem];
        }
        domain.m_fx[gnode] += fx;
        domain.m_fy[gnode] += fy;
        domain.m_fz[gnode] += fz;
    }

    Release(&fz_elem);
	Release(&fy_elem);
	Release(&fx_elem);

    printf("fx: ");
    for( Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode )
        printf("%f, ", scalar<Real_t>(domain.m_fx[gnode]));
    printf("\n");
    printf("fy: ");
    for( Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode )
        printf("%f, ", scalar<Real_t>(domain.m_fy[gnode]));
    printf("\n");
    printf("fz: ");
    for( Index_t gnode=0 ; gnode<domain.numNode() ; ++gnode )
        printf("%f, ", scalar<Real_t>(domain.m_fz[gnode]));
    printf("\n");
}

static inline
void VoluDer(const Real_t x0, const Real_t x1, const Real_t x2,
             const Real_t x3, const Real_t x4, const Real_t x5,
             const Real_t y0, const Real_t y1, const Real_t y2,
             const Real_t y3, const Real_t y4, const Real_t y5,
             const Real_t z0, const Real_t z1, const Real_t z2,
             const Real_t z3, const Real_t z4, const Real_t z5,
             Real_t* dvdx, Real_t* dvdy, Real_t* dvdz) {
	const Real_t twelfth = Real_t(1.0) / Real_t(12.0) ;
    *dvdx =
	(y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +
	(y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -
	(y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5);
	*dvdy =
	- (x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) -
	(x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +
	(x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5);

	*dvdz =
	- (y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) -
	(y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +
	(y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5);

	*dvdx *= twelfth;
	*dvdy *= twelfth;
	*dvdz *= twelfth;
}

static inline
void CalcElemVolumeDerivative(Real_t dvdx[8],
                              Real_t dvdy[8],
                              Real_t dvdz[8],
                              const Real_t x[8],
                              const Real_t y[8],
                              const Real_t z[8]) {
	VoluDer(x[1], x[2], x[3], x[4], x[5], x[7],
			y[1], y[2], y[3], y[4], y[5], y[7],
			z[1], z[2], z[3], z[4], z[5], z[7],
			&dvdx[0], &dvdy[0], &dvdz[0]);
	VoluDer(x[0], x[1], x[2], x[7], x[4], x[6],
			y[0], y[1], y[2], y[7], y[4], y[6],
			z[0], z[1], z[2], z[7], z[4], z[6],
			&dvdx[3], &dvdy[3], &dvdz[3]);
	VoluDer(x[3], x[0], x[1], x[6], x[7], x[5],
			y[3], y[0], y[1], y[6], y[7], y[5],
			z[3], z[0], z[1], z[6], z[7], z[5],
			&dvdx[2], &dvdy[2], &dvdz[2]);
	VoluDer(x[2], x[3], x[0], x[5], x[6], x[4],
			y[2], y[3], y[0], y[5], y[6], y[4],
			z[2], z[3], z[0], z[5], z[6], z[4],
			&dvdx[1], &dvdy[1], &dvdz[1]);
	VoluDer(x[7], x[6], x[5], x[0], x[3], x[1],
			y[7], y[6], y[5], y[0], y[3], y[1],
			z[7], z[6], z[5], z[0], z[3], z[1],
			&dvdx[4], &dvdy[4], &dvdz[4]);
	VoluDer(x[4], x[7], x[6], x[1], x[0], x[2],
			y[4], y[7], y[6], y[1], y[0], y[2],
			z[4], z[7], z[6], z[1], z[0], z[2],
			&dvdx[5], &dvdy[5], &dvdz[5]);
	VoluDer(x[5], x[4], x[7], x[2], x[1], x[3],
			y[5], y[4], y[7], y[2], y[1], y[3],
			z[5], z[4], z[7], z[2], z[1], z[3],
			&dvdx[6], &dvdy[6], &dvdz[6]);
	VoluDer(x[6], x[5], x[4], x[3], x[2], x[0],
			y[6], y[5], y[4], y[3], y[2], y[0],
			z[6], z[5], z[4], z[3], z[2], z[0],
			&dvdx[7], &dvdy[7], &dvdz[7]);
}

static inline
void CalcHourglassControlForElems(multi_array<Real_t>& determ, Real_t hgcoef) {
	Index_t numElem = domain.numElem();
    Real_t *dvdx = Allocate<Real_t>(numElem * 8);
	Real_t *dvdy = Allocate<Real_t>(numElem * 8);
	Real_t *dvdz = Allocate<Real_t>(numElem * 8);
	Real_t *x8n  = Allocate<Real_t>(numElem * 8);
	Real_t *y8n  = Allocate<Real_t>(numElem * 8);
	Real_t *z8n  = Allocate<Real_t>(numElem * 8);

//TODO: Optimize

    for (Index_t i=0 ; i<numElem ; ++i) {
        Real_t  x1[8],  y1[8],  z1[8];
        Real_t pfx[8], pfy[8], pfz[8];

        CollectDomainNodesToElemNodes(i, x1, y1, z1);

        CalcElemVolumeDerivative(pfx, pfy, pfz, x1, y1, z1);

		for(Index_t ii=0;ii<8;++ii) {
			Index_t jj=8*i+ii;

			dvdx[jj] = pfx[ii];
			dvdy[jj] = pfy[ii];
			dvdz[jj] = pfz[ii];
			x8n[jj]  = x1[ii];
			y8n[jj]  = y1[ii];
			z8n[jj]  = z1[ii];
		}

		determ[i] = domain.m_volo[i] * domain.m_v[i];
		if (scalar<Real_t>(domain.m_v[i]) <= Real_t(0.0))
			exit(VolumeError) ;
    }

    if (hgcoef > Real_t(0.))
		CalcFBHourglassForceForElems(determ,x8n,y8n,z8n,dvdx,dvdy,dvdz,hgcoef);

    Release(&z8n);
	Release(&y8n);
	Release(&x8n);
	Release(&dvdz);
	Release(&dvdy);
	Release(&dvdx);
}

static inline
void CalcVolumeForceForElems() {
    Index_t numElem = domain.numElem();
    if (numElem != 0) {
        multi_array<Real_t> sigxx;
        multi_array<Real_t> sigyy;
        multi_array<Real_t> sigzz;
        multi_array<Real_t> determ;

        sigxx = zeros<Real_t>(numElem);
        sigyy = zeros<Real_t>(numElem);
        sigzz = zeros<Real_t>(numElem);
        determ = zeros<Real_t>(numElem);

        printf("InitStressTermsForElems\n");
        InitStressTermsForElems(sigxx, sigyy, sigzz);

        printf("IntegrateStressForElems\n");
        IntegrateStressForElems(sigxx, sigyy, sigzz,
                                determ, domain.numElem());

        if (scalar(sum(as<Real_t>(determ <= Real_t(0.0)))))
            exit(VolumeError);

        printf("CalcHourglassControlForElems\n");
        CalcHourglassControlForElems(determ, domain.hgcoef());
    }
}

static inline
void CalcForceForNodes() {
    Index_t numNode = domain.numNode();

    domain.m_fx = zeros<Real_t>(numNode);
    domain.m_fy = zeros<Real_t>(numNode);
    domain.m_fz = zeros<Real_t>(numNode);

    /* Calcforce calls partial, force, hourq */
    printf("CalcVolumeForceForElems\n");
    CalcVolumeForceForElems();
}

static inline
void LagrangeNodal() {
    const Real_t delt = domain.deltatime() ;

    printf("CalcForceForNodes\n");
    CalcForceForNodes();

    printf("CalcAccelerationForNodes\n");
    CalcAccelerationForNodes();

    printf("ApplyAccelerationBoundaryConditionsForNodes\n");
    ApplyAccelerationBoundaryConditionsForNodes();

    printf("CalcVelocityForNodes\n");
    CalcVelocityForNodes(delt, domain.u_cut());

    printf("CalcPositionForNodes\n");
    CalcPositionForNodes(delt);

    return;
}

static inline
Real_t AreaFace(const Real_t x0, const Real_t x1,
				const Real_t x2, const Real_t x3,
				const Real_t y0, const Real_t y1,
				const Real_t y2, const Real_t y3,
				const Real_t z0, const Real_t z1,
				const Real_t z2, const Real_t z3) {
	Real_t fx = (x2 - x0) - (x3 - x1);
	Real_t fy = (y2 - y0) - (y3 - y1);
	Real_t fz = (z2 - z0) - (z3 - z1);
	Real_t gx = (x2 - x0) + (x3 - x1);
	Real_t gy = (y2 - y0) + (y3 - y1);
	Real_t gz = (z2 - z0) + (z3 - z1);
	Real_t area =
	(fx * fx + fy * fy + fz * fz) *
	(gx * gx + gy * gy + gz * gz) -
	(fx * gx + fy * gy + fz * gz) *
	(fx * gx + fy * gy + fz * gz);
	return area ;
}

static inline
Real_t CalcElemCharacteristicLength(const Real_t x[8],
									const Real_t y[8],
									const Real_t z[8],
									const Real_t volume) {
	Real_t a, charLength = Real_t(0.0);

	a = AreaFace(x[0],x[1],x[2],x[3],
				 y[0],y[1],y[2],y[3],
				 z[0],z[1],z[2],z[3]) ;
	charLength = std::max(a,charLength) ;

	a = AreaFace(x[4],x[5],x[6],x[7],
				 y[4],y[5],y[6],y[7],
				 z[4],z[5],z[6],z[7]) ;
	charLength = std::max(a,charLength) ;

	a = AreaFace(x[0],x[1],x[5],x[4],
				 y[0],y[1],y[5],y[4],
				 z[0],z[1],z[5],z[4]) ;
	charLength = std::max(a,charLength) ;

	a = AreaFace(x[1],x[2],x[6],x[5],
				 y[1],y[2],y[6],y[5],
				 z[1],z[2],z[6],z[5]) ;
	charLength = std::max(a,charLength) ;

	a = AreaFace(x[2],x[3],x[7],x[6],
				 y[2],y[3],y[7],y[6],
				 z[2],z[3],z[7],z[6]) ;
	charLength = std::max(a,charLength) ;

	a = AreaFace(x[3],x[0],x[4],x[7],
				 y[3],y[0],y[4],y[7],
				 z[3],z[0],z[4],z[7]) ;
	charLength = std::max(a,charLength) ;

	charLength = Real_t(4.0) * volume / SQRT(charLength);

	return charLength;
}

static inline
void CalcElemVelocityGrandient(const Real_t* const xvel,
							   const Real_t* const yvel,
							   const Real_t* const zvel,
							   const Real_t b[][8],
							   const Real_t detJ,
							   Real_t* const d) {
	const Real_t inv_detJ = Real_t(1.0) / detJ ;
	Real_t dyddx, dxddy, dzddx, dxddz, dzddy, dyddz;
	const Real_t* const pfx = b[0];
	const Real_t* const pfy = b[1];
	const Real_t* const pfz = b[2];

	d[0] = inv_detJ * ( pfx[0] * (xvel[0]-xvel[6])
					   + pfx[1] * (xvel[1]-xvel[7])
					   + pfx[2] * (xvel[2]-xvel[4])
					   + pfx[3] * (xvel[3]-xvel[5]) );

	d[1] = inv_detJ * ( pfy[0] * (yvel[0]-yvel[6])
					   + pfy[1] * (yvel[1]-yvel[7])
					   + pfy[2] * (yvel[2]-yvel[4])
					   + pfy[3] * (yvel[3]-yvel[5]) );

	d[2] = inv_detJ * ( pfz[0] * (zvel[0]-zvel[6])
					   + pfz[1] * (zvel[1]-zvel[7])
					   + pfz[2] * (zvel[2]-zvel[4])
					   + pfz[3] * (zvel[3]-zvel[5]) );

	dyddx  = inv_detJ * ( pfx[0] * (yvel[0]-yvel[6])
						 + pfx[1] * (yvel[1]-yvel[7])
						 + pfx[2] * (yvel[2]-yvel[4])
						 + pfx[3] * (yvel[3]-yvel[5]) );

	dxddy  = inv_detJ * ( pfy[0] * (xvel[0]-xvel[6])
						 + pfy[1] * (xvel[1]-xvel[7])
						 + pfy[2] * (xvel[2]-xvel[4])
						 + pfy[3] * (xvel[3]-xvel[5]) );

	dzddx  = inv_detJ * ( pfx[0] * (zvel[0]-zvel[6])
						 + pfx[1] * (zvel[1]-zvel[7])
						 + pfx[2] * (zvel[2]-zvel[4])
						 + pfx[3] * (zvel[3]-zvel[5]) );

	dxddz  = inv_detJ * ( pfz[0] * (xvel[0]-xvel[6])
						 + pfz[1] * (xvel[1]-xvel[7])
						 + pfz[2] * (xvel[2]-xvel[4])
						 + pfz[3] * (xvel[3]-xvel[5]) );

	dzddy  = inv_detJ * ( pfy[0] * (zvel[0]-zvel[6])
						 + pfy[1] * (zvel[1]-zvel[7])
						 + pfy[2] * (zvel[2]-zvel[4])
						 + pfy[3] * (zvel[3]-zvel[5]) );

	dyddz  = inv_detJ * ( pfz[0] * (yvel[0]-yvel[6])
						 + pfz[1] * (yvel[1]-yvel[7])
						 + pfz[2] * (yvel[2]-yvel[4])
						 + pfz[3] * (yvel[3]-yvel[5]) );
	// d[5]  = Real_t( .5) * ( dxddy + dyddx );
	// d[4]  = Real_t( .5) * ( dxddz + dzddx );
	// d[3]  = Real_t( .5) * ( dzddy + dyddz );
}

static inline
void CalcKinematicsForElems(Index_t numElem, Real_t dt) {
//TODO: Optimize
    for(Index_t k=0 ; k<numElem ; ++k) {
        Real_t B[3][8];
        Real_t D[6];
		Real_t x_local[8];
		Real_t y_local[8];
		Real_t z_local[8];
		Real_t xd_local[8];
		Real_t yd_local[8];
		Real_t zd_local[8];
		Real_t detJ = Real_t(0.0);

        Real_t volume, relativeVolume;

        for( Index_t lnode=0 ; lnode<8 ; ++lnode ) {
			Index_t gnode = scalar<Index_t>(domain.m_nodelist[Index_t(8)*k+lnode]);
			x_local[lnode] = scalar<Real_t>(domain.m_x[gnode]);
			y_local[lnode] = scalar<Real_t>(domain.m_y[gnode]);
			z_local[lnode] = scalar<Real_t>(domain.m_z[gnode]);
		}

        volume = CalcElemVolume(x_local, y_local, z_local);
        relativeVolume = volume / scalar<Real_t>(domain.m_volo[k]);
        domain.m_vnew[k] = relativeVolume;
		domain.m_delv[k] = relativeVolume - domain.m_v[k];

        domain.m_arealg[k] = CalcElemCharacteristicLength(x_local, y_local, z_local, volume);

        for( Index_t lnode=0 ; lnode<8 ; ++lnode ) {
			Index_t gnode = scalar<Index_t>(domain.m_nodelist[Index_t(8)*k+lnode]);
			xd_local[lnode] = scalar<Real_t>(domain.m_xd[gnode]);
			yd_local[lnode] = scalar<Real_t>(domain.m_yd[gnode]);
			zd_local[lnode] = scalar<Real_t>(domain.m_zd[gnode]);
		}

        Real_t dt2 = Real_t(0.5) * dt;
		for ( Index_t j=0 ; j<8 ; ++j ) {
			x_local[j] -= dt2 * xd_local[j];
			y_local[j] -= dt2 * yd_local[j];
			z_local[j] -= dt2 * zd_local[j];
		}

        detJ = CalcElemShapeFunctionDerivatives(x_local, y_local, z_local, B);
        // printf("detJ: %d\n", detJ);

		CalcElemVelocityGrandient(xd_local, yd_local, zd_local, B, detJ, D);

        // printf("D: \t[0] = %f, \t[1] = %f, \t[2] = %f\n", D[0],D[1],D[2]);
        domain.m_dxx[k] = D[0];
		domain.m_dyy[k] = D[1];
		domain.m_dzz[k] = D[2];
    }
}

static inline
void CalcLagrangeElements(Real_t deltatime) {
    Index_t numElem = domain.numElem() ;
    if (numElem > 0) {
        CalcKinematicsForElems(numElem, deltatime);

        multi_array<Real_t> vdov;
        vdov = domain.m_dxx+domain.m_dyy+domain.m_dzz;
        // std::cout << "dxx: " <<  domain.m_dxx << std::endl;
        // std::cout << "dyy: " <<  domain.m_dyy << std::endl;
        // std::cout << "dzz: " <<  domain.m_dzz << std::endl;

        multi_array<Real_t> vdovthird;
        vdovthird = vdov/Real_t(3.0);

        domain.m_vdov(vdov);
        // std::cout << "vdov: " <<  domain.m_vdov << std::endl;
        domain.m_dxx -= vdovthird;
        domain.m_dyy -= vdovthird;
        domain.m_dzz -= vdovthird;

        if (scalar(sum(as<Real_t>(domain.m_vnew <= Real_t(0.0)))))
            exit(VolumeError);
    }
}

static inline
void CalcMonotonicQGradientsForElems() {
//TODO: Optimize
#define SUM4(a,b,c,d) (a + b + c + d)
	for (Index_t i = 0 ; i < domain.numElem() ; ++i ) {
		const Real_t ptiny = Real_t(1.e-36) ;
		Real_t ax,ay,az;
		Real_t dxv,dyv,dzv;

		Index_t n0 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+0]);
		Index_t n1 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+1]);
		Index_t n2 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+2]);
		Index_t n3 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+3]);
		Index_t n4 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+4]);
		Index_t n5 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+5]);
		Index_t n6 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+6]);
		Index_t n7 = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+7]);

		Real_t x0 = scalar<Real_t>(domain.m_x[n0]);
		Real_t x1 = scalar<Real_t>(domain.m_x[n1]);
		Real_t x2 = scalar<Real_t>(domain.m_x[n2]);
		Real_t x3 = scalar<Real_t>(domain.m_x[n3]);
		Real_t x4 = scalar<Real_t>(domain.m_x[n4]);
		Real_t x5 = scalar<Real_t>(domain.m_x[n5]);
		Real_t x6 = scalar<Real_t>(domain.m_x[n6]);
		Real_t x7 = scalar<Real_t>(domain.m_x[n7]);

		Real_t y0 = scalar<Real_t>(domain.m_y[n0]);
		Real_t y1 = scalar<Real_t>(domain.m_y[n1]);
		Real_t y2 = scalar<Real_t>(domain.m_y[n2]);
		Real_t y3 = scalar<Real_t>(domain.m_y[n3]);
		Real_t y4 = scalar<Real_t>(domain.m_y[n4]);
		Real_t y5 = scalar<Real_t>(domain.m_y[n5]);
		Real_t y6 = scalar<Real_t>(domain.m_y[n6]);
		Real_t y7 = scalar<Real_t>(domain.m_y[n7]);

		Real_t z0 = scalar<Real_t>(domain.m_z[n0]);
		Real_t z1 = scalar<Real_t>(domain.m_z[n1]);
		Real_t z2 = scalar<Real_t>(domain.m_z[n2]);
		Real_t z3 = scalar<Real_t>(domain.m_z[n3]);
		Real_t z4 = scalar<Real_t>(domain.m_z[n4]);
		Real_t z5 = scalar<Real_t>(domain.m_z[n5]);
		Real_t z6 = scalar<Real_t>(domain.m_z[n6]);
		Real_t z7 = scalar<Real_t>(domain.m_z[n7]);

		Real_t xv0 = scalar<Real_t>(domain.m_xd[n0]);
		Real_t xv1 = scalar<Real_t>(domain.m_xd[n1]);
		Real_t xv2 = scalar<Real_t>(domain.m_xd[n2]);
		Real_t xv3 = scalar<Real_t>(domain.m_xd[n3]);
		Real_t xv4 = scalar<Real_t>(domain.m_xd[n4]);
		Real_t xv5 = scalar<Real_t>(domain.m_xd[n5]);
		Real_t xv6 = scalar<Real_t>(domain.m_xd[n6]);
		Real_t xv7 = scalar<Real_t>(domain.m_xd[n7]);

		Real_t yv0 = scalar<Real_t>(domain.m_yd[n0]);
		Real_t yv1 = scalar<Real_t>(domain.m_yd[n1]);
		Real_t yv2 = scalar<Real_t>(domain.m_yd[n2]);
		Real_t yv3 = scalar<Real_t>(domain.m_yd[n3]);
		Real_t yv4 = scalar<Real_t>(domain.m_yd[n4]);
		Real_t yv5 = scalar<Real_t>(domain.m_yd[n5]);
		Real_t yv6 = scalar<Real_t>(domain.m_yd[n6]);
		Real_t yv7 = scalar<Real_t>(domain.m_yd[n7]);

		Real_t zv0 = scalar<Real_t>(domain.m_zd[n0]);
		Real_t zv1 = scalar<Real_t>(domain.m_zd[n1]);
		Real_t zv2 = scalar<Real_t>(domain.m_zd[n2]);
		Real_t zv3 = scalar<Real_t>(domain.m_zd[n3]);
		Real_t zv4 = scalar<Real_t>(domain.m_zd[n4]);
		Real_t zv5 = scalar<Real_t>(domain.m_zd[n5]);
		Real_t zv6 = scalar<Real_t>(domain.m_zd[n6]);
		Real_t zv7 = scalar<Real_t>(domain.m_zd[n7]);

		Real_t vol = scalar<Real_t>(domain.m_volo[i]*domain.m_vnew[i]);
		Real_t norm = Real_t(1.0) / ( vol + ptiny );

		Real_t dxj = Real_t(-0.25)*(SUM4(x0,x1,x5,x4) - SUM4(x3,x2,x6,x7)) ;
		Real_t dyj = Real_t(-0.25)*(SUM4(y0,y1,y5,y4) - SUM4(y3,y2,y6,y7)) ;
		Real_t dzj = Real_t(-0.25)*(SUM4(z0,z1,z5,z4) - SUM4(z3,z2,z6,z7)) ;

		Real_t dxi = Real_t( 0.25)*(SUM4(x1,x2,x6,x5) - SUM4(x0,x3,x7,x4)) ;
		Real_t dyi = Real_t( 0.25)*(SUM4(y1,y2,y6,y5) - SUM4(y0,y3,y7,y4)) ;
		Real_t dzi = Real_t( 0.25)*(SUM4(z1,z2,z6,z5) - SUM4(z0,z3,z7,z4)) ;

		Real_t dxk = Real_t( 0.25)*(SUM4(x4,x5,x6,x7) - SUM4(x0,x1,x2,x3)) ;
		Real_t dyk = Real_t( 0.25)*(SUM4(y4,y5,y6,y7) - SUM4(y0,y1,y2,y3)) ;
		Real_t dzk = Real_t( 0.25)*(SUM4(z4,z5,z6,z7) - SUM4(z0,z1,z2,z3)) ;

		/* find delvk and delxk ( i cross j ) */

		ax = dyi*dzj - dzi*dyj ;
		ay = dzi*dxj - dxi*dzj ;
		az = dxi*dyj - dyi*dxj ;

		domain.m_delx_zeta[i] = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

		ax *= norm ;
		ay *= norm ;
		az *= norm ;

		dxv = Real_t(0.25)*(SUM4(xv4,xv5,xv6,xv7) - SUM4(xv0,xv1,xv2,xv3)) ;
		dyv = Real_t(0.25)*(SUM4(yv4,yv5,yv6,yv7) - SUM4(yv0,yv1,yv2,yv3)) ;
		dzv = Real_t(0.25)*(SUM4(zv4,zv5,zv6,zv7) - SUM4(zv0,zv1,zv2,zv3)) ;

		domain.m_delv_zeta[i] = ax*dxv + ay*dyv + az*dzv ;

		/* find delxi and delvi ( j cross k ) */

		ax = dyj*dzk - dzj*dyk ;
		ay = dzj*dxk - dxj*dzk ;
		az = dxj*dyk - dyj*dxk ;

		domain.m_delx_xi[i] = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

		ax *= norm ;
		ay *= norm ;
		az *= norm ;

		dxv = Real_t(0.25)*(SUM4(xv1,xv2,xv6,xv5) - SUM4(xv0,xv3,xv7,xv4)) ;
		dyv = Real_t(0.25)*(SUM4(yv1,yv2,yv6,yv5) - SUM4(yv0,yv3,yv7,yv4)) ;
		dzv = Real_t(0.25)*(SUM4(zv1,zv2,zv6,zv5) - SUM4(zv0,zv3,zv7,zv4)) ;

		domain.m_delv_xi[i] = ax*dxv + ay*dyv + az*dzv ;

		/* find delxj and delvj ( k cross i ) */

		ax = dyk*dzi - dzk*dyi ;
		ay = dzk*dxi - dxk*dzi ;
		az = dxk*dyi - dyk*dxi ;

		domain.m_delx_eta[i] = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

		ax *= norm ;
		ay *= norm ;
		az *= norm ;

		dxv = Real_t(-0.25)*(SUM4(xv0,xv1,xv5,xv4) - SUM4(xv3,xv2,xv6,xv7)) ;
		dyv = Real_t(-0.25)*(SUM4(yv0,yv1,yv5,yv4) - SUM4(yv3,yv2,yv6,yv7)) ;
		dzv = Real_t(-0.25)*(SUM4(zv0,zv1,zv5,zv4) - SUM4(zv3,zv2,zv6,zv7)) ;

		domain.m_delv_eta[i] = ax*dxv + ay*dyv + az*dzv ;
	}
#undef SUM4
}

static inline
void CalcMonotonicQRegionForElems(Real_t qlc_monoq,
                                  Real_t qqc_monoq,
                                  Real_t monoq_limiter_mult,
                                  Real_t monoq_max_slope,
                                  Real_t ptiny,
                                  Index_t elength) {
    /* phixi */
    multi_array<Real_t> norm;
    norm = Real_t(1.)/(domain.m_delv_xi+ptiny);

    multi_array<Real_t> delvm, delvp;
    delvm = zeros<Real_t>(elength);
    delvp = zeros<Real_t>(elength);

    multi_array<Int_t> bcMask;
    bcMask = domain.m_elemBC & XI_M;
    delvm(as<Real_t>(bcMask == 0)*gather(delvm, domain.m_delv_xi, domain.m_lxim)
                + as<Real_t>(bcMask == XI_M_SYMM)*domain.m_delv_xi
                + as<Real_t>(bcMask == XI_M_FREE)*Real_t(0.0));

    bcMask(domain.m_elemBC & XI_P);
    delvp(as<Real_t>(bcMask == 0)*gather(delvp, domain.m_delv_xi, domain.m_lxip)
                + as<Real_t>(bcMask == XI_P_SYMM)*domain.m_delv_xi
                + as<Real_t>(bcMask == XI_P_FREE)*Real_t(0.0));
    // delvm[_(1, -1)] = as<Real_t>(bcMask == 0 || bcMask == XI_M_SYMM)
    //                     * domain.m_delv_xi[_(0, -2)];
    // delvm[0] = as<Real_t>(bcMask == 0 || bcMask == XI_M_SYMM)*domain.m_delv_xi[0];
    //
    // bcMask = domain.m_elemBC & XI_P;
    // delvp[_(0, -2)] = as<Real_t>(bcMask == 0 || bcMask == XI_P_SYMM)
    //                     * domain.m_delv_xi[_(0, -2)];
    // delvp[-1] = as<Real_t>(bcMask == 0 || bcMask == XI_P_SYMM)
    //                 * domain.m_delv_xi[-1];

    delvm(delvm*norm);
    delvp(delvp*norm);

    multi_array<Real_t> phixi;
    phixi = Real_t(.5)*(delvm+delvp);

    delvm *= monoq_limiter_mult;
    delvp *= monoq_limiter_mult;

    phixi(as<Real_t>(delvm < phixi)*delvm + as<Real_t>(delvm >= phixi)*phixi);
    phixi(as<Real_t>(delvp < phixi)*delvp + as<Real_t>(delvp >= phixi)*phixi);
    phixi(as<Real_t>(phixi >= Real_t(0.))*phixi);
    phixi(as<Real_t>(phixi > monoq_max_slope)*monoq_max_slope + as<Real_t>(phixi <= monoq_max_slope)*phixi);

    /* phieta */
    // Real_t sizeX = domain.sizeX();
    norm(Real_t(1.)/(domain.m_delv_eta+ptiny));
    bcMask(domain.m_elemBC & ETA_M);
    delvm(as<Real_t>(bcMask == 0)*gather(delvm, domain.m_delv_eta, domain.m_letam)
                + as<Real_t>(bcMask == ETA_M_SYMM)*domain.m_delv_eta
                + as<Real_t>(bcMask == ETA_M_FREE)*Real_t(0.0));

    bcMask(domain.m_elemBC & ETA_P);
    delvp(as<Real_t>(bcMask == 0)*gather(delvp, domain.m_delv_eta, domain.m_letap)
                + as<Real_t>(bcMask == ETA_P_SYMM)*domain.m_delv_eta
                + as<Real_t>(bcMask == ETA_P_FREE)*Real_t(0.0));
    // delvm[_(0, sizeX, 1)] = as<Real_t>(bcMask == 0)*domain.m_delv_eta[_(0, sizeX, 1)]
    //             + as<Real_t>(bcMask == ETA_M_SYMM)*domain.m_delv_eta
    //             + as<Real_t>(bcMask == ETA_M_FREE)*Real_t(0.0);
    // delvm[_(sizeX, -1, 1)] = as<Real_t>(bcMask == 0)*domain.m_delv_eta[_(0, domain.m_delv_eta.len()-sizeX, 1)]
    //             + as<Real_t>(bcMask == ETA_M_SYMM)*domain.m_delv_eta
    //             + as<Real_t>(bcMask == ETA_M_FREE)*Real_t(0.0);
    //
    // bcMask = domain.m_elemBC & ETA_P;
    // delvp = as<Real_t>(bcMask == 0)*domain.m_delv_eta[_(domain.sizeX(), -1, 1)]
    //             + as<Real_t>(bcMask == ETA_P_SYMM)*domain.m_delv_eta
    //             + as<Real_t>(bcMask == ETA_P_FREE)*Real_t(0.0);

    delvm(delvm*norm);
    delvp(delvp*norm);

    multi_array<Real_t> phieta;
    phieta = Real_t(.5)*(delvm+delvp);

    delvm *= monoq_limiter_mult;
    delvp *= monoq_limiter_mult;

    phieta(as<Real_t>(delvm < phieta)*delvm + as<Real_t>(delvm >= phieta)*phieta);
    phieta(as<Real_t>(delvp < phieta)*delvp + as<Real_t>(delvp >= phieta)*phieta);
    phieta(as<Real_t>(phieta >= Real_t(0.))*phieta);
    phieta(as<Real_t>(phieta > monoq_max_slope)*monoq_max_slope + as<Real_t>(phieta <= monoq_max_slope)*phieta);

    /*  phizeta */
    // Real_t sizeX2 = sizeX*sizeX;
    norm(Real_t(1.)/(domain.m_delv_zeta+ptiny));
    bcMask(domain.m_elemBC & ZETA_M);
    delvm(as<Real_t>(bcMask == 0)*gather(delvm, domain.m_delv_zeta, domain.m_lzetam)
                + as<Real_t>(bcMask == ZETA_M_SYMM)*domain.m_delv_zeta
                + as<Real_t>(bcMask == ZETA_M_FREE)*Real_t(0.0));

    bcMask(domain.m_elemBC & ZETA_P);
    delvp(as<Real_t>(bcMask == 0)*gather(delvp, domain.m_delv_zeta, domain.m_lzetap)
                + as<Real_t>(bcMask == ZETA_P_SYMM)*domain.m_delv_zeta
                + as<Real_t>(bcMask == ZETA_P_FREE)*Real_t(0.0));

    // bcMask = domain.m_elemBC & ZETA_M;
    // delvm[_(0, sizeX2, 1)] = as<Real_t>(bcMask == 0)*domain.m_delv_zeta[_(0, sizeX2, 1)]
    //             + as<Real_t>(bcMask == ZETA_M_SYMM)*domain.m_delv_zeta
    //             + as<Real_t>(bcMask == ZETA_M_FREE)*Real_t(0.0);
    // delvm[_(sizeX2, -1, 1)] = as<Real_t>(bcMask == 0)*domain.m_delv_zeta[_(0, domain.m_delv_zeta.len()-sizeX2, 1)]
    //             + as<Real_t>(bcMask == ZETA_M_SYMM)*domain.m_delv_zeta
    //             + as<Real_t>(bcMask == ZETA_M_FREE)*Real_t(0.0);
    //
    // bcMask = domain.m_elemBC & ZETA_P;
    // Real_t slice = domain.numElem() - sizeX2;
    // delvp[_(0, slice, 1)] = as<Real_t>(bcMask == 0)*domain.m_delv_zeta[_(sizeX2, slice+sizeX2, 1)]
    //             + as<Real_t>(bcMask == ZETA_P_SYMM)*domain.m_delv_zeta
    //             + as<Real_t>(bcMask == ZETA_P_FREE)*Real_t(0.0);
    // delvp[_(slice, -1, 1)] = as<Real_t>(bcMask == 0)*domain.m_delv_zeta[_(slice, -1, 1)]
    //             + as<Real_t>(bcMask == ZETA_P_SYMM)*domain.m_delv_zeta
    //             + as<Real_t>(bcMask == ZETA_P_FREE)*Real_t(0.0);

    delvm(delvm*norm);
    delvp(delvp*norm);

    multi_array<Real_t> phizeta;
    phizeta = Real_t(.5)*(delvm+delvp);

    delvm *= monoq_limiter_mult;
    delvp *= monoq_limiter_mult;

    phizeta(as<Real_t>(delvm < phizeta)*delvm + as<Real_t>(delvm >= phizeta)*phizeta);
    phizeta(as<Real_t>(delvp < phizeta)*delvp + as<Real_t>(delvp >= phizeta)*phizeta);
    phizeta(as<Real_t>(phizeta >= Real_t(0.))*phizeta);
    phizeta(as<Real_t>(phizeta > monoq_max_slope)*monoq_max_slope + as<Real_t>(phizeta <= monoq_max_slope)*phizeta);

    /* Remove length scale */
    multi_array<Real_t> delvxxi;
    multi_array<Real_t> delvxeta;
    multi_array<Real_t> delvxzeta;
    delvxxi = domain.m_delv_xi*domain.m_delx_xi;
    delvxxi(as<Real_t>(delvxxi <= Real_t(0.))*delvxxi);

    delvxeta = domain.m_delv_eta*domain.m_delx_eta;
    delvxeta(as<Real_t>(delvxeta <= Real_t(0.))*delvxeta);

    delvxzeta = domain.m_delv_zeta*domain.m_delx_zeta;
    delvxzeta(as<Real_t>(delvxzeta <= Real_t(0.))*delvxzeta);

    multi_array<Real_t> rho;
    rho = domain.m_elemMass/(domain.m_volo*domain.m_vnew);

    multi_array<Real_t> qlin;
    qlin = -qlc_monoq*rho*(delvxxi*(Real_t(1.)-phixi)
        +delvxeta*(Real_t(1.)-phieta)+delvxzeta*(Real_t(1.)-phizeta));
    domain.m_ql(as<Real_t>(domain.m_vdov > Real_t(0.))*Real_t(0.)
                    + as<Real_t>(domain.m_vdov <= Real_t(0.))*qlin);

    multi_array<Real_t> qquad;
    qquad = qqc_monoq*rho*(delvxxi*delvxxi*(Real_t(1.)-phixi*phixi)
        +delvxeta*delvxeta*(Real_t(1.)-phieta*phieta)
        +delvxzeta*delvxzeta*(Real_t(1.)-phizeta*phizeta));
    domain.m_qq(as<Real_t>(domain.m_vdov > Real_t(0.))*Real_t(0.)
                    + as<Real_t>(domain.m_vdov <= Real_t(0.))*qquad);
}

static inline
void CalcMonotonicQForElems() {
    if (domain.numElem() > 0)
        CalcMonotonicQRegionForElems(domain.qlc_monoq(),
                                     domain.qqc_monoq(),
                                     domain.monoq_limiter_mult(),
                                     domain.monoq_max_slope(),
                                     Real_t(1.e-36),
                                     domain.numElem());
}

static inline
void CalcQForElems() {
    /* Calculate velocity gradients */
    printf("CalcMonotonicQGradientsForElems\n");
    CalcMonotonicQGradientsForElems() ;

    /* Transfer veloctiy gradients in the first order elements */
    printf("CalcMonotonicQForElems\n");
    CalcMonotonicQForElems() ;

    /* Don't allow excessive artificial viscosity */
    Index_t numElem = domain.numElem();

    if (numElem != 0) {
        if (scalar(sum(as<Real_t>(domain.m_q > domain.qstop()))))
            exit(QStopError);
    }
}

static inline
void CalcSoundSpeedForElems(multi_array<Real_t>& vnewc, Real_t rho0,
                            multi_array<Real_t>& enewc, multi_array<Real_t>& pnewc,
                            multi_array<Real_t>& pbvc, multi_array<Real_t>& bvc) {
    multi_array<Real_t> ssTmp;
    ssTmp = (pbvc * enewc + vnewc * vnewc * bvc * pnewc) / rho0;
    ssTmp(as<Real_t>(domain.m_ss <= Real_t(.111111e-36))*Real_t(.333333e-18)
                + as<Real_t>(domain.m_ss > Real_t(.111111e-36))*sqrt(ssTmp));
    domain.m_ss(ssTmp);
}

static inline
void CalcPressureForElems(multi_array<Real_t>& p_new, multi_array<Real_t>& bvc,
                          multi_array<Real_t>& pbvc, multi_array<Real_t>& e_old,
                          multi_array<Real_t>& compression, multi_array<Real_t>& vnewc,
                          Real_t pmin, Real_t p_cut, Real_t eosvmax) {
    Real_t cls = Real_t(2.0)/Real_t(3.0);
    bvc = cls*(compression+Real_t(1.));
    pbvc(cls);
    p_new(bvc*e_old);
    p_new(as<Real_t>(abs(p_new) >= p_cut)*p_new);
    p_new(as<Real_t>(vnewc < eosvmax)*p_new);
    p_new(as<Real_t>(p_new >= pmin)*p_new + as<Real_t>(p_new < pmin)*pmin);
}

static inline
void CalcEnergyForElems(multi_array<Real_t>& p_new, multi_array<Real_t>& e_new,
                        multi_array<Real_t>& q_new, multi_array<Real_t>& bvc,
                        multi_array<Real_t>& pbvc, multi_array<Real_t>& p_old,
                        multi_array<Real_t>& e_old, multi_array<Real_t>& q_old,
                        multi_array<Real_t>& compression, multi_array<Real_t>& compHalfStep,
                        multi_array<Real_t>& vnewc, multi_array<Real_t>& work,
                        multi_array<Real_t>& delvc, Real_t pmin,
                        Real_t p_cut, Real_t  e_cut, Real_t q_cut, Real_t emin,
                        multi_array<Real_t>& qq, multi_array<Real_t>& ql,
                        Real_t rho0,
                        Real_t eosvmax) {
    multi_array<Real_t> pHalfStep;
    pHalfStep = zeros<Real_t>(domain.numElem());

    // printf("CalcEnergyForElems\n");
    e_new(e_old-Real_t(0.5)*delvc*(p_old+q_old)+Real_t(0.5)*work);
    e_new(as<Real_t>(e_new < emin)*emin + as<Real_t>(e_new >= emin)*e_new);

    CalcPressureForElems(pHalfStep, bvc, pbvc, e_new, compHalfStep,
                         vnewc, pmin, p_cut, eosvmax);

    multi_array<Real_t> vhalf;
    vhalf = Real_t(1.) / (Real_t(1.) + compHalfStep);
    multi_array<Real_t> ssc;
    ssc = (pbvc*e_new+vhalf*vhalf*bvc*pHalfStep)/rho0;
    ssc(as<Real_t>(ssc <= Real_t(.111111e-36))*Real_t(.333333e-18)
            + as<Real_t>(ssc > Real_t(.111111e-36))*sqrt(ssc));

    q_new(as<Real_t>(delvc <= Real_t(0.))*(ssc*ql+qq));

    e_new(e_new+Real_t(0.5)*delvc * (Real_t(3.0)*(p_old+q_old)
                -Real_t(4.0)*(pHalfStep+q_new)));

    e_new += Real_t(0.5)*work;
    e_new(as<Real_t>(abs(e_new) >= e_cut)*e_new);
    e_new(as<Real_t>(e_new < emin)*emin + as<Real_t>(e_new >= emin)*e_new);

    CalcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc,
						 pmin, p_cut, eosvmax);

    ssc((pbvc*e_new+vnewc*vnewc*bvc*p_new)/rho0);
    ssc(as<Real_t>(ssc <= Real_t(.111111e-36))*Real_t(.333333e-18)
            + as<Real_t>(ssc > Real_t(.111111e-36))*sqrt(ssc));

    multi_array<Real_t> q_tilde;
    q_tilde = as<Real_t>(delvc > Real_t(0.))*Real_t(0.)
                + as<Real_t>(delvc <= Real_t(0.))*(ssc*ql+qq);

    e_new(e_new-(Real_t(7.0)*(p_old+q_old)-Real_t(8.0)*(pHalfStep+q_new)
                +(p_new+q_tilde))*delvc*(Real_t(1.0)/Real_t(6.0)));
    e_new(as<Real_t>(abs(e_new) >= e_cut)*e_new);
    e_new(as<Real_t>(e_new < emin)*emin + as<Real_t>(e_new >= emin)*e_new);

    CalcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc,
						 pmin, p_cut, eosvmax);

    if (scalar(sum(as<Real_t>(delvc <= Real_t(0.))))) {
//TODO: Maybe wrong
        ssc((pbvc*e_new+vnewc*vnewc*bvc*p_new)/rho0);
        ssc(as<Real_t>(ssc <= Real_t(.111111e-36))*Real_t(.333333e-18)
                + as<Real_t>(ssc > Real_t(.111111e-36))*sqrt(ssc));

        q_new(ssc*ql+qq);
        q_new(as<Real_t>(abs(q_new >= q_cut))*q_new);
    }
}

static inline
void EvalEOSForElems(multi_array<Real_t>& vnewc, Index_t length) {
    // printf("EvalEOSForElems\n");
    multi_array<Real_t> e_old;
    e_old = domain.m_e;

    multi_array<Real_t> delvc;
    delvc = domain.m_delv;

    multi_array<Real_t> p_old;
    p_old = domain.m_p;

    multi_array<Real_t> q_old;
    q_old = domain.m_q;

    multi_array<Real_t> compression;
    multi_array<Real_t> compHalfStep;
    compression = Real_t(1.) / vnewc - Real_t(1.);
    compHalfStep = Real_t(1.) / (vnewc - delvc * Real_t(.5)) - Real_t(1.);

    if (domain.eosvmin() != Real_t(0.))
        compHalfStep(as<Real_t>(vnewc <= domain.eosvmin())*compression
                        + as<Real_t>(vnewc > domain.eosvmin())*compHalfStep);

    if (domain.eosvmax() != Real_t(0.)) {
        p_old = as<Real_t>(vnewc < domain.eosvmax())*p_old;
        compression(as<Real_t>(vnewc < domain.eosvmax())*compression);
        compHalfStep(as<Real_t>(vnewc < domain.eosvmax())*compHalfStep);
    }

    multi_array<Real_t> work;
    work = zeros<Real_t>(length);

    multi_array<Real_t> qq;
    qq = domain.m_qq;

    multi_array<Real_t> ql;
    ql = domain.m_ql;


    multi_array<Real_t> p_new;
    p_new = zeros<Real_t>(length);
    multi_array<Real_t> e_new;
    e_new = zeros<Real_t>(length);
    multi_array<Real_t> q_new;
    q_new = zeros<Real_t>(length);
    multi_array<Real_t> bvc;
    bvc = zeros<Real_t>(length);
    multi_array<Real_t> pbvc;
    pbvc = zeros<Real_t>(length);

    printf("CalcEnergyForElems\n");
    CalcEnergyForElems(p_new, e_new, q_new, bvc, pbvc,
  					   p_old, e_old,  q_old, compression, compHalfStep,
  					   vnewc, work,  delvc, domain.pmin(),
                       domain.p_cut(), domain.e_cut(), domain.q_cut(), domain.emin(),
  					   qq, ql, domain.refdens(), domain.eosvmax());

    domain.m_p(p_new);
    domain.m_e(e_new);
    domain.m_q(q_new);

    printf("CalcSoundSpeedForElems\n");
    CalcSoundSpeedForElems(vnewc, domain.refdens(), e_new, p_new, pbvc, bvc);
}

static inline
void ApplyMaterialPropertiesForElems() {
    Index_t length = domain.numElem() ;

    if (length != 0) {
        Real_t eosvmin = domain.eosvmin();
        Real_t eosvmax = domain.eosvmax();

        multi_array<Real_t> vnewc;
        vnewc = domain.m_vnew;

        if (eosvmin != Real_t(0.)) {
            vnewc = as<Real_t>(vnewc < eosvmin)*eosvmin +
                      as<Real_t>(vnewc >= eosvmin)*vnewc;
            domain.m_v(as<Real_t>(domain.m_v < eosvmin)*eosvmin
                          + as<Real_t>(domain.m_v >= eosvmin)*domain.m_v);

        }
        if (eosvmax != Real_t(0.)) {
            vnewc = as<Real_t>(vnewc > eosvmax)*eosvmax +
                      as<Real_t>(vnewc <= eosvmax)*vnewc;
            domain.m_v(as<Real_t>(domain.m_v > eosvmax)*eosvmax
                          + as<Real_t>(domain.m_v <= eosvmax)*domain.m_v);
        }

        if (scalar(sum(as<Real_t>(domain.m_v <= 0.))))
            exit(VolumeError);

        printf("EvalEOSForElems\n");
        EvalEOSForElems(vnewc, length);
    }
}

static inline
void UpdateVolumesForElems() {
    if (domain.numElem() != 0)
        domain.m_v(as<Real_t>(abs(domain.m_vnew - Real_t(1.0)) >= domain.v_cut())*domain.m_vnew);
}

static inline
void LagrangeElements() {
    printf("CalcLagrangeElements\n");
    CalcLagrangeElements(domain.deltatime());

    printf("CalcQForElems\n");
    CalcQForElems();

    printf("ApplyMaterialPropertiesForElems\n");
    ApplyMaterialPropertiesForElems();

    printf("UpdateVolumesForElems\n");
    UpdateVolumesForElems();
}

static inline
void CalcCourantConstraintForElems()
{
    Real_t dtcourant = Real_t(1.0e+20) ;
    Index_t   courant_elem = -1 ;
    Real_t      qqc = domain.qqc() ;
    Index_t length = domain.numElem() ;

    Real_t  qqc2 = Real_t(64.0) * qqc * qqc ;

    for (Index_t i = 0 ; i < length ; ++i) {
        Index_t indx = i;

        Real_t dtf = scalar<Real_t>(domain.m_ss[indx] * domain.m_ss[indx]);

        if (scalar<bool>(domain.m_vdov[indx] < Real_t(0.))) {

            dtf = scalar<Real_t>(dtf
            + qqc2 * domain.m_arealg[indx] * domain.m_arealg[indx]
            * domain.m_vdov[indx] * domain.m_vdov[indx]);
        }

        dtf = SQRT(dtf) ;

        dtf = scalar<Real_t>(domain.m_arealg[indx] / dtf);

        /* determine minimum timestep with its corresponding elem */
        if (scalar<bool>(domain.m_vdov[indx] != Real_t(0.))) {
            if ( dtf < dtcourant ) {
                dtcourant = dtf ;
                courant_elem = indx ;
            }
        }
    }

    /* Don't try to register a time constraint if none of the elements
     * were active */
    if (courant_elem != -1) {
        domain.dtcourant() = dtcourant ;
        printf("dtcourant: %f\n", domain.dtcourant());
    }

    return ;
}

static inline
void CalcHydroConstraintForElems()
{
    Real_t dthydro = Real_t(1.0e+20) ;
    Index_t hydro_elem = -1 ;
    Real_t dvovmax = domain.dvovmax() ;
    Index_t length = domain.numElem() ;

    for (Index_t i = 0 ; i < length ; ++i) {
        Index_t indx = i;

        if (scalar<bool>(domain.m_vdov[indx] != Real_t(0.))) {
            Real_t dtdvov = dvovmax / (FABS(scalar<Real_t>(domain.m_vdov[indx])+Real_t(1.e-20))) ;
            if ( dthydro > dtdvov ) {
                dthydro = dtdvov ;
                hydro_elem = indx ;
            }
        }
    }

    if (hydro_elem != -1) {
        domain.dthydro() = dthydro ;
        printf("dthydro: %f\n\n", domain.dthydro());
    }

    return ;
}

// static inline
// void CalcCourantConstraintForElems() {
//     Real_t qqc = domain.qqc() ;
// 	Real_t qqc2 = Real_t(64.0) * qqc * qqc ;
//
//     multi_array<Real_t> dtf;
//     dtf = domain.m_ss*domain.m_ss;
//
//     multi_array<Real_t> dtfConst;
//     dtfConst = dtf+qqc2*domain.m_arealg*domain.m_arealg*domain.m_vdov*domain.m_vdov;
//     dtf = as<Real_t>(domain.m_vdov < Real_t(0.))*dtfConst
//             + as<Real_t>(domain.m_vdov >= Real_t(0.))*dtf;
//
//     dtf = sqrt(dtf);
//     dtf = domain.m_arealg/dtf;
//     dtf = as<Real_t>(domain.m_vdov != Real_t(0.))*dtf
//             + as<Real_t>(domain.m_vdov == Real_t(0.))*Real_t(1.e+20);
//
//     Real_t minDtf = scalar<Real_t>(min(dtf));
//
//     if (minDtf < Real_t(1.0e+20)) {
//         domain.dtcourant() = minDtf;
//         std::cout << "Values: " <<  dtf << std::endl;
//         std::cout << "Minimum: " << bh_minimum_reduce(dtf, (int64_t)0) << std::endl;
//         printf("dtcourant: %f\n", domain.dtcourant());
//     }
// }
//
// static inline
// void CalcHydroConstraintForElems() {
//     if (scalar(as<Real_t>(domain.m_vdov != Real_t(0.)))) {
//         multi_array<Real_t> dtdvov;
//         dtdvov = domain.dvovmax()/(abs(domain.m_vdov)+Real_t(1.e-20));
//         dtdvov = as<Real_t>(domain.m_vdov != Real_t(0.))*dtdvov
//                     + as<Real_t>(domain.m_vdov == Real_t(0.))*Real_t(1.e+20);
//         Real_t minDtdvov = scalar<Real_t>(min(dtdvov));
//
//         if (minDtdvov < Real_t(1.0e+20)) {
//             domain.dthydro() = minDtdvov;
//
//             std::cout << "Values: " <<  dtdvov << std::endl;
//             std::cout << "Minimum: " << bh_minimum_reduce(dtdvov, (int64_t)-1) << std::endl;
//             printf("dthydro: %f\n", domain.dthydro());
//         }
//     }
// }

static inline
void CalcTimeConstraintsForElems() {
    /* evaluate time constraint */
    CalcCourantConstraintForElems() ;

    /* check hydro constraint */
    CalcHydroConstraintForElems() ;
}

static inline
void LagrangeLeapFrog() {
    /* calculate nodal forces, accelerations, velocities, positions, with
     * applied boundary conditions and slide surface considerations */
    LagrangeNodal();
    printf("LagrangeElements\n");

    /* calculate element quantities (i.e. velocity gradient & q), and update
     * material states */
    LagrangeElements();

    CalcTimeConstraintsForElems();
}

static inline
void TimeIncrement() {
    printf("C: %f, \tH: %f\n", domain.dtcourant(), domain.dthydro());

    if ((domain.dtfixed() <= Real_t(0.0))
            && (domain.cycle() != Int_t(0))) {
        Real_t ratio;

        /* This will require a reduction in parallel */
        Real_t newdt = Real_t(1.0e+20);
        if (domain.dtcourant() < newdt)
            newdt = domain.dtcourant() / Real_t(2.0);

        if (domain.dthydro() < newdt)
            newdt = domain.dthydro() * Real_t(2.0) / Real_t(3.0);

        Real_t olddt = domain.deltatime();
        ratio = newdt / olddt;
        if (ratio >= Real_t(1.0)) {
            if (ratio < domain.deltatimemultlb())
                newdt = olddt;
            else if (ratio > domain.deltatimemultub())
                newdt = olddt*domain.deltatimemultub();
        }

        if (newdt > domain.dtmax())
            newdt = domain.dtmax();
        domain.deltatime() = newdt ;
    }

    Real_t targetdt = domain.stoptime() - domain.time();

    /* TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE */
    if ((targetdt > domain.deltatime()) &&
            (targetdt < (Real_t(4.0) * domain.deltatime() / Real_t(3.0))))
        targetdt = Real_t(2.0) * domain.deltatime() / Real_t(3.0) ;

    if (targetdt < domain.deltatime())
        domain.deltatime() = targetdt;

    domain.time() += domain.deltatime();

    printf("time = %e, dt = %e, stop = %e\n",
           double(domain.time()), double(domain.deltatime()), double(domain.stoptime()));
    ++domain.cycle();
}

int main(int argc, char *argv[]){
    printf("Main\n");
    Index_t edgeElems = atoi(argv[1]);
    Index_t edgeNodes = edgeElems+1;

    Real_t tx, ty, tz;
    Index_t nidx, zidx;

    /* Initialize Sedov Mesh */
    domain.sizeX() = edgeElems;
	domain.sizeY() = edgeElems;
	domain.sizeZ() = edgeElems;
	domain.numElem() = edgeElems*edgeElems*edgeElems;
	domain.numNode() = edgeNodes*edgeNodes*edgeNodes;

    Index_t domElems = domain.numElem();

    /* allocate field memory */
    domain.AllocateElemPersistent(domElems);
    domain.AllocateElemTemporary(domElems);
    domain.AllocateNodalPersistent(domain.numNode());
    domain.AllocateNodesets(edgeNodes*edgeNodes);

    /* initialize nodal coordinates */
    nidx = 0 ;
	tz  = Real_t(0.) ;
	for (Index_t plane=0; plane<edgeNodes; ++plane) {
		ty = Real_t(0.) ;
		for (Index_t row=0; row<edgeNodes; ++row) {
			tx = Real_t(0.) ;
			for (Index_t col=0; col<edgeNodes; ++col) {
				domain.m_x[nidx] = tx;
				domain.m_y[nidx] = ty;
				domain.m_z[nidx] = tz;
				 ++nidx ;
				tx = Real_t(1.125)*Real_t(col+1)/Real_t(edgeElems) ;
			}
			ty = Real_t(1.125)*Real_t(row+1)/Real_t(edgeElems) ;
		}
		tz = Real_t(1.125)*Real_t(plane+1)/Real_t(edgeElems) ;

	}
    /* embed hexehedral elements in nodal point lattice */
    nidx = 0;
    zidx = 0;
    for (Index_t plane=0; plane<edgeElems; ++plane) {
        for (Index_t row=0; row<edgeElems; ++row) {
            for (Index_t col=0; col<edgeElems; ++col) {
                Index_t pos = Index_t(8)*zidx;
                domain.m_nodelist[pos](nidx);
                domain.m_nodelist[pos+1] = nidx+1;
                domain.m_nodelist[pos+2] = nidx+edgeNodes+1;
                domain.m_nodelist[pos+3] = nidx+edgeNodes;
                domain.m_nodelist[pos+4] = nidx+edgeNodes*edgeNodes;
                domain.m_nodelist[pos+5] = nidx+edgeNodes*edgeNodes+1;
                domain.m_nodelist[pos+6] = nidx+edgeNodes*edgeNodes+edgeNodes+1;
                domain.m_nodelist[pos+7] = nidx+edgeNodes*edgeNodes+edgeNodes;
                ++zidx;
                ++nidx;
            }
            ++nidx;
        }
        nidx += edgeNodes;
    }

    domain.AllocateNodeElemIndexes();

    /* initialize material parameters */
    domain.dtfixed() = Real_t(-1.0e-7);
    domain.deltatime() = Real_t(1.0e-7);
    domain.deltatimemultlb() = Real_t(1.1);
    domain.deltatimemultub() = Real_t(1.2);
    domain.stoptime() = Real_t(1.0e-2);
    domain.dtcourant() = Real_t(1.0e+20);
    domain.dthydro() = Real_t(1.0e+20);
    domain.dtmax() = Real_t(1.0e-2);
    domain.time() = Real_t(0.);
    domain.cycle() = 0;

    domain.e_cut() = Real_t(1.0e-7);
    domain.p_cut() = Real_t(1.0e-7);
    domain.q_cut() = Real_t(1.0e-7);
    domain.u_cut() = Real_t(1.0e-7);
    domain.v_cut() = Real_t(1.0e-10);

    domain.hgcoef() = Real_t(3.0);
    domain.ss4o3() = Real_t(4.0)/Real_t(3.0);

    domain.qstop() = Real_t(1.0e+12);
    domain.monoq_max_slope() = Real_t(1.0);
    domain.monoq_limiter_mult() = Real_t(2.0);
    domain.qlc_monoq() = Real_t(0.5);
    domain.qqc_monoq() = Real_t(2.0)/Real_t(3.0);
    domain.qqc() = Real_t(2.0);

    domain.pmin() = Real_t(0.);
    domain.emin() = Real_t(-1.0e+15);

    domain.dvovmax() = Real_t(0.1);

    domain.eosvmax() = Real_t(1.0e+9);
    domain.eosvmin() = Real_t(1.0e-9);

    domain.refdens() = Real_t(1.0);

    /* initialize field data */
    for (Index_t i=0; i<domElems; ++i) {
        Real_t x_local[8], y_local[8], z_local[8];
        for(Index_t lnode=0; lnode<8; ++lnode) {
            Index_t gnode = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+lnode]);
            x_local[lnode] = scalar<Real_t>(domain.m_x[gnode]);
            y_local[lnode] = scalar<Real_t>(domain.m_y[gnode]);
            z_local[lnode] = scalar<Real_t>(domain.m_z[gnode]);
        }

        printf("x_local: ");
        for(Index_t lnode=0; lnode<8; ++lnode)
            printf("%f, ", x_local[lnode]);
        printf("\n");
        printf("y_local: ");
        for(Index_t lnode=0; lnode<8; ++lnode)
            printf("%f, ", y_local[lnode]);
        printf("\n");
        printf("z_local: ");
        for(Index_t lnode=0; lnode<8; ++lnode)
            printf("%f, ", z_local[lnode]);
        printf("\n");

        // volume calculations
        Real_t volume = CalcElemVolume(x_local, y_local, z_local);
        printf("Volumen: %d\n", volume);
        domain.m_volo[i] = volume;
        domain.m_elemMass[i] = volume;
        for (Index_t j=0; j<8; ++j) {
            Index_t idx = scalar<Index_t>(domain.m_nodelist[Index_t(8)*i+j]);
            domain.m_nodalMass[idx] += volume / Real_t(8.0);
        }
    }

    /* deposit energy */
    domain.m_e[0] = Real_t(3.948746e+7);

    nidx = 0 ;
	for (Index_t i=0; i<edgeNodes; ++i) {
		Index_t planeInc = i*edgeNodes*edgeNodes;
		Index_t rowInc   = i*edgeNodes;
		for (Index_t j=0; j<edgeNodes; ++j) {
			domain.m_symmX[nidx] = planeInc + j*edgeNodes;
			domain.m_symmY[nidx] = planeInc + j;
			domain.m_symmZ[nidx] = rowInc + j;
			++nidx ;
		}
	}

    domain.m_lxim[0] = 0;
	for (Index_t i=1; i<domElems; ++i) {
		domain.m_lxim[i] = i-1 ;
		domain.m_lxip[i-1] = i ;
	}
	domain.m_lxip[domElems-1] = domElems-1 ;

    for (Index_t i=0; i<edgeElems; ++i) {
        domain.m_letam[i] = i;
        domain.m_letap[domElems-edgeElems+i] = domElems-edgeElems+i ;
    }

    for (Index_t i=edgeElems; i<domElems; ++i) {
        domain.m_letam[i] = i-edgeElems ;
        domain.m_letap[i-edgeElems] = i ;
    }

    for (Index_t i=0; i<edgeElems*edgeElems; ++i) {
        domain.m_lzetam[i] = i;
        domain.m_lzetap[domElems-edgeElems*edgeElems+i] = domElems-edgeElems*edgeElems+i ;
    }

    for (Index_t i=edgeElems*edgeElems; i<domElems; ++i) {
        domain.m_lzetam[i] = i - edgeElems*edgeElems ;
        domain.m_lzetap[i-edgeElems*edgeElems] = i ;
    }

    /* set up boundary condition information */
    domain.m_elemBC = zeros<Int_t>(domElems);

    /* faces on "external" boundaries will be */
    /* symmetry plane or free surface BCs */
    for (Index_t i=0; i<edgeElems; ++i) {
        Index_t planeInc = i*edgeElems*edgeElems ;
        Index_t rowInc = i*edgeElems ;
        for (Index_t j=0; j<edgeElems; ++j) {
            domain.m_elemBC[planeInc+j*edgeElems] |= XI_M_SYMM ;
            domain.m_elemBC[planeInc+j*edgeElems+edgeElems-1] |= XI_P_FREE ;
            domain.m_elemBC[planeInc+j] |= ETA_M_SYMM ;
            domain.m_elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_FREE ;
            domain.m_elemBC[rowInc+j] |= ZETA_M_SYMM ;
            domain.m_elemBC[rowInc+j+domElems-edgeElems*edgeElems] |= ZETA_P_FREE ;
        }
    }

    timeval start, end;
	gettimeofday(&start, NULL);

    printf("Starting\n");
    while(domain.time() < domain.stoptime() ) {
        if (domain.cycle() == 306)
			break;
        TimeIncrement();
        // printf("LagrangeLeapFrog\n");
        LagrangeLeapFrog();

        printf("--------------\n\n Wave: %d\n", domain.cycle());

#if LULESH_SHOW_PROGRESS
        printf("time = %e, dt=%e\n",
               double(domain.time()), double(domain.deltatime()) ) ;
#endif
    }

    gettimeofday(&end, NULL);
	double elapsed_time = double(end.tv_sec - start.tv_sec) + double(end.tv_usec - start.tv_usec) *1e-6;

	printf("\n\nElapsed time = %12.6e\n\n", elapsed_time);

	Index_t ElemId = 0;
	printf("Run completed:\n");
	printf("   Problem size        =  %i \n", edgeElems);
	printf("   Iteration count     =  %i \n", domain.cycle());
	printf("   Final Origin Energy = %12.6e \n", scalar<Real_t>(domain.m_e[ElemId]));

	Real_t   MaxAbsDiff = Real_t(0.0);
	Real_t TotalAbsDiff = Real_t(0.0);
	Real_t   MaxRelDiff = Real_t(0.0);

	for (Index_t j=0; j<edgeElems; ++j) {
        for (Index_t k=j+1; k<edgeElems; ++k) {
            Real_t AbsDiff = FABS(scalar<Real_t>(domain.m_e[j*edgeElems+k])
                    - scalar<Real_t>(domain.m_e[k*edgeElems+j]));
            TotalAbsDiff += AbsDiff;

            if (MaxAbsDiff < AbsDiff)
                MaxAbsDiff = AbsDiff;

            Real_t RelDiff = AbsDiff / scalar<Real_t>(domain.m_e[k*edgeElems+j]);
            if (MaxRelDiff < RelDiff)
                MaxRelDiff = RelDiff;
        }
    }

    printf("Testing Plane 0 of Energy Array:\n");
    printf("   MaxAbsDiff   = %12.6e\n", MaxAbsDiff);
    printf("   TotalAbsDiff = %12.6e\n", TotalAbsDiff);
    printf("   MaxRelDiff   = %12.6e\n\n", MaxRelDiff);
	return 0;
}
