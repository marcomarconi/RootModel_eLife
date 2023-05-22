#ifndef PBD_HPP
#define PBD_HPP
#include <Process.hpp>
#include "tissue.hpp"
#include <Attributes.hpp>
#include <Function.hpp>
#include <MDXProcessTissue.hpp>
#include <MeshProcessSystem.hpp>
#include <Contour.hpp>
#include <CCTopology.hpp>
#include <Mesh.hpp>
#include <MeshBuilder.hpp>
#include <Process.hpp>
#include <MDXProcessCellDivide.hpp>
#include <MeshProcessStructure.hpp>
#include <ToolsProcessSystem.hpp>
#include <CCVerify.hpp>
#include <CCIndex.hpp>
#include <CellMakerMesh2D.hpp>
#include <GeomMathUtils.hpp>
#include <Geometry.hpp>
#include <Matrix.hpp>
#include <Triangulate.hpp>

using namespace mdx;

class PBD : public Process{
public:
    PBD(const Process& process)
        : Process(process) {
        setName("Model/Root/07 PBD Engine");

        addParm("Print Stats",
                "Print Stats",
                "True",
                QStringList() << "True"
                              << "False");
        addParm("Substrate Fixed",
                "Substrate Fixed",
                "True",
                QStringList() << "True"
                              << "False");
        addParm("Velocity Update",
                "Solver Order",
                "First Order",
                QStringList() << "First Order"
                              << "Second Order");
        addParm("Max Velocity", "", "5");
        addParm("Velocity Static Damping", "", "1");
        addParm("Viscous Damping", "", "10");
        addParm("Momentum Preservation K", "", "0");

        addParm("PBD Paramters", "PBD Paramters", "");
        addParm("PBD iterations", "PBD iterations", "5");
        addParm("PBD stiffness correction", "PBD stiffness correction",  "True", QStringList() << "False"
                              << "True");
        addParm("PBD pressure stiffness", "stiffness area", "1");
        addParm("PBD pressure alpha", "PBD pressure alpha", "10");
        addParm("PBD distance stiffness", "stiffness distance", "1");
        addParm("PBD distance alpha", "PBD distance alpha", "0");
        addParm("PBD bending stiffness", "stiffness bending", "0");
        addParm("PBD strain stiffness", "PBD stiffness area", "1");
        //addParm("PBD strain alpha", "PBD strain alpha", "0");
        addParm("PBD strain shear stiffness", "PBD strain shear stiffness", "1");
        addParm("PBD strain normalize", "PBD strain normalize", "True",
                                                                        QStringList() << "False"
                                                                                      << "True");
        addParm("PBD shape stiffness", "PBD stiffness shape", "1");
        addParm("PBD shape allow stretching", "PBD shape allow stretchinge", "True",
                                                                        QStringList() << "False"
                                                                                      << "True");
    }



    static bool init_ShapeMatchingConstraint(const std::vector<Point3d> x0, const std::vector<double> invMasses,
                                      int numPoints,
                                      Point3d& restCm,
                                      Matrix3d& invRestMat) ;

    bool solve_ShapeMatchingConstraint(const std::vector<Point3d> x0,
                                       const std::vector<Point3d> x,
                                       const std::vector<double> invMasses,
                                       int numPoints,
                                       const Point3d& restCm,
                                       const Matrix3d& invRestMat,
                                       bool allowStretch,
                                       std::vector<Point3d>& corr,
                                       Matrix3d* rot) ;
    bool solve_DistanceConstraint(
        const Point3d &p0, double invMass0,
        const Point3d &p1, double invMass1,
        const double restLength,
        const double compressionStiffness,
        const double stretchStiffness,
        Point3d &corr0, Point3d &corr1);


    bool solve_DistanceConstraint_XPBD(const Point3d &p0, double invMass0,
                                   const Point3d &p1, double invMass1,
                                   const double restLength,
                                   const double compressionStiffness,
                                   const double stretchStiffness,
                                   Point3d &corr0, Point3d &corr1, double &lambda, double compliance) ;

    // ----------------------------------------------------------------------------------------------
    static bool init_StrainTriangleConstraint(
        const Point3d &p0,
        const Point3d &p1,
        const Point3d &p2,
        Matrix2d &invRestMat);

    // ----------------------------------------------------------------------------------------------
    bool solve_StrainTriangleConstraint(
            const Point3d &p0,
            const Point3d &p1,
            const Point3d &p2,
            const Matrix2d &invRestMat,
            const double xxStiffness,
            const double yyStiffness,
            const double xyStiffness,
            const bool normalizeStretch,
            const bool normalizeShear,
            Point3d &corr0, Point3d &corr1, Point3d &corr2, Point3d &lambda, double alpha_tilde);

    bool solve_StrainTriangleConstraint(
             const Point3d &p0, double invMass0,
             const Point3d &p1, double invMass1,
             const Point3d &p2, double invMass2,
             const Matrix2d &invRestMat,
             const double xxStiffness,
             const double yyStiffness,
             const double xyStiffness,
             const bool normalizeStretch,
             const bool normalizeShear,
             Point3d &corr0, Point3d &corr1, Point3d &corr2);

    bool solve_PressureConstraint(
           const CCStructure& cs,
           const CCIndexDataAttr &indexAttr,
           Tissue::VertexDataAttr &vMAttr,
           const std::set<CCIndex> fs,
           const double restArea,
           const double pressure,
           std::map<CCIndex, Point3d>& corr, double &lambda, double compliance);

    bool solve_DihedralConstraint(
         const Point3d &p0, double invMass0,
         const Point3d &p1, double invMass1,
         const Point3d &p2, double invMass2,
         const Point3d &p3, double invMass3,
         const double restAngle,
         const double stiffness,
         Point3d &corr0, Point3d &corr1, Point3d &corr2, Point3d &corr3, bool verbose);

    bool init_FEMTriangleConstraint(
            const Point3d &p0,
            const Point3d &p1,
            const Point3d &p2,
            double &area,
            Matrix2d &invRestMat);

    bool solve_FEMTriangleConstraint(const Point3d &p0, double invMass0,
            const Point3d &p1, double invMass1,
            const Point3d &p2, double invMass2,
            const double &area,
            const Matrix2d &invRestMat,
            const double youngsModulusX,
            const double youngsModulusY,
            const double youngsModulusShear,
            const double poissonRatioXY,
            const double poissonRatioYX,
            Point3d &corr0, Point3d &corr1, Point3d &corr2);



    void momentumPreservation(const CCStructure& cs, const CCIndexDataAttr &indexAttr, Tissue::VertexDataAttr &vMAttr, double K);
    void dampVelocity(CCIndex v, const CCStructure& cs, const CCIndexDataAttr& indexAttr,
                      Tissue::VertexDataAttr &vMAttr, Tissue::VertexData& vD );
    void semiImplicitEuler(CCIndex v, const CCStructure& cs,
                            CCIndexDataAttr& indexAttr, Tissue::VertexDataAttr &vMAttr, Tissue::VertexData& vD);
    void solve();

    bool initialize(QWidget* parent);
    bool rewind(QWidget* parent);
    void update(double Dt);

    bool printStats = false;
    double convergenceLagParm = 0;
    QString velocityUpdate;
    double maxVelocity = 5;
    double staticDamping = 0, viscousDamping = 0;
    int PBDiterations = 0;
    bool stiffnessCorrection = true;
    double PBDpressureStiffness = 0;
    double PBDpressureK = 0;
    double PBDpressureAlpha = 0;
    double PBDdistanceStiffness = 0;
    double PBDbendingStiffness = 0;
    double PBDdistanceAlpha = 0;
    double PBDstrainStiffness = 0;
    double PBDstrainAlpha = 0;
    double PBDstrainShearStiffness = 0;
    bool PBDstrainNormalize = true;
    double PBDshapeStiffness = 0;
    bool PBDshapeAllowStretching = true;

private:
    double Dt = 0;
    CCIndexDataAttr *indexAttr = 0;
    int debug_step = 0;

};


#endif // PBD_HPP
