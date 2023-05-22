//
// This file is part of MorphoDynamX - http://www.MorphoDynamX.org
// Copyright (C) 2012-2016 Richard S. Smith and collaborators.
//
// If you use MorphoDynamX in your work, please cite:
//   http://dx.doi.org/10.7554/eLife.05864
//
// MorphoDynamX is free software, and is licensed under under the terms of the
// GNU General (GPL) Public License version 2.0, http://www.gnu.org/licenses.
// beeeeeeeeeeeeeeeh

#ifndef ROOT_HPP
#define ROOT_HPP

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "PBD.hpp"
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
#include <CellMakerUtils.hpp>

#include <CCVerify.hpp>
#include <CCIndex.hpp>
#include <CellMakerMesh2D.hpp>
#include <GeomMathUtils.hpp>
#include <Geometry.hpp>
#include <Matrix.hpp>
#include <Triangulate.hpp>
#include <limits>
#include <math.h>
#include <CCDivideCell.hpp>
#include <cstdlib>
#include "tbb/concurrent_unordered_map.h"
#include <fstream>
#include <chrono>
#include <complex>
#include <ctime>

using namespace mdx;

namespace ROOT
{

// pre-definition
class Root;
class Remesh;
class MechanicalGrowth;
class SetGlobalAttr;
class ClearCells;
class TriangulateFacesX;
class DeleteCell;
class ExecuteTest;
class AddFace;




class Mechanics : public Process{
public:
    Mechanics(const Process& process)
        : Process(process) {
        setName("Model/Root/05 Mechanics");
        addParm("Dt", "Time step in hours", "0.03");
        addParm("Converge Threshold", "Threshold for convergence", ".1");
        addParm("Convergence Lag", "Convergence Lag", "1");
        addParm("Verbose", "Verbose", "True",
                QStringList() << "True"
                              << "False");
        addParm("Debug Steps", "Debug Steps", "10");

        // Mass Springs
        addParm("Mass Spring Parameters", "", "");
        addParm("Wall EK", "Stiffness of the cross springs", "1");
        addParm("Wall CK", "Stiffness of the cross springs", "1");
        addParm("Shear EK", "Extensional Stiffness of the shear edge springs", "0");
        addParm("Shear CK", "Compression Stiffness of the shear edge springs", "1");
        addParm("Auxin-induced wall relaxation K1", "Auxin-induced wall relaxation K1", "0.05");
        addParm("Auxin-induced wall relaxation K2", "Auxin-induced wall relaxation K2", "3");
        addParm("Wall stress", "Wall stress", "1");
        addParm("Wall stress K1", "Wall stress K1", "0.01");
        addParm("Wall stress K2", "Wall stress K2", "2");
        // Hydrostatics
        addParm("Hydrostatic Parameters", "", "");
        addParm("Turgor Pressure", "Value of the turgor pressure in the cells", "2");
        addParm("Turgor Pressure Rate", "Value of the rate of turgor pressure in the cells", "0.5");
        // External Forces
        addParm("External Forces", "", "");
        addParm("Gravity Force", "Gravity Force", "0");
        addParm("Gravity Direction", "Gravity Direction", "0,-1,0");
        addParm("Friction", "Friction", "0");
        // Misc
        addParm("Tissue Process", "Name of Tissue Process", "Model/Root/03 Cell Tissue");
        addParm("PBD Engine",
                "Name of PBD Engine",
                "Model/Root/07 PBD Engine");
    }


    double calcNorm();
    bool initialize(QWidget* parent);
    bool rewind(QWidget* parent);
    bool step();
    PBD* PBDProcess = 0;
    double Dt = 0;
    double userTime = 0;
    double realTime = 0;
    std::vector<double> growthRatesVector;

private:
    void calcForces(CCIndex v, const CCStructure& cs, const CCIndexDataAttr& indexAttr,
                    Tissue::CellDataAttr &cellAttr, Tissue::FaceDataAttr &faceAttr,
                    Tissue::EdgeDataAttr &edgeAttr, Tissue::VertexData& vD);

    Point3d calcForcesFace(CCIndex f,
        Tissue::CellData& cD, Tissue::FaceData& fD, Tissue::EdgeData& eD, Tissue::VertexData& vD);
    Point3d calcForcesEdge(
        const CCStructure& cs, const CCIndexDataAttr& indexAttr, CCIndex e, CCIndex v, int label);

    CCIndexDataAttr* indexAttr = 0;
    Tissue* tissueProcess = 0;
    double wallStress = 0, wallStressK1 = 0, wallStressK2 = 0;
    Point3d gravity;
    double convergeThresh = 1e-6;
    double convergenceLag = 0;
    double normal = 0, prev_normal;
    Point3d prevQCcentroid , QCcentroid;
    int debug_step = 0;
    bool Verbose = false;
};

class SetDirichlet : public Process {
public:
    SetDirichlet(const Process& process)
        : Process(process) {

        setName("Model/Root/24 Set Dirichlet");
        setDesc("Assign Dirichlet to vertexes.");
        setIcon(QIcon(":/images/CellType.png"));

        addParm("Dirichlet", "Set Dirichlet condition on selected vtx in x,y,z", "0 0 0");
    }
    bool initialize(QWidget* parent);
    bool run();

private:
    Point3u dirichletFixed;
};


class MechanicalGrowth : public Process {
public:
    MechanicalGrowth(const Process& process)
        : Process(process) {
        setName("Model/Root/052 Mechanical Growth");
        addParm("Verbose", "Verbose", "True",
                QStringList() << "True"
                              << "False");
        addParm("Polarity Method",
                "Polarity Method",
                "Cell Axis",
                QStringList() << "Cell Axis"
                              << "Principal Strains"
                              << "Membrane Stretch");
        addParm("Strain Tensor",
                "Strain Tensor",
                "Green Strain Tensor",
                QStringList() << "Green Strain Tensor"
                              << "Shape Strain Tensor");
        addParm("Walls Growth", "Walls Growth", "");
        addParm("Strain Threshold for Growth", "Strain Threshold for Growth", "0.01");
        addParm("Walls Growth Rate", "Walls Growth Rate", "10");
        addParm("Area Growth Rate", "Area Growth Rate", "0");
        addParm("MicroFibrils", "MicroFibrils", "");
        addParm("MF reorientation rate", "MF Reorientation Rate", "0.02");
        addParm("MF reorientation strainrate",
                "MF Reorientation strainrate",
                "0.02");
        addParm("MF Degradation", "MF Degradation", "0.01");
        addParm("MF Delete After Division", "MF Delete After Division","True",
                QStringList() << "True"
                              << "False");
        addParm("Zonation", "Zonation", "");
        addParm("Elongation Zone", "Elongation Zone", "50");
        addParm("Differentiation Zone", "Differentiation Zone", "100");
        addParm("Mechanics Process",
                "Name of Mechanics derivatives process",
                "Model/Root/05 Mechanics");
        addParm("Mechanical Solver Process",
                "Name of Mechanical Solver Process",
                "Model/Root/051 Mechanical Solver");

    }

    bool initialize(QWidget* parent);
    bool step(double Dt);

private:
    Mechanics* mechanicsProcess = 0;
    PBD* PBDProcess = 0;
    bool Verbose = false;
};


/* Suggestions:
 * - High auxin -> higher auxin degradation
 * - High auxin -> lower PIN degradation on membranes  (Paciorek et al. 2005.
 * - High auxin -> Lower PIN expression  (Vieten 2005) or higher PIN degradation in cytoplasm (Baster et al., 2013).
 * - Growth dynamically affect by auxin? (tested: does not seems to work well)
*/
class Chemicals : public Process {
public:
    Chemicals(const Process& process)
        : Process(process) {
        setName("Model/Root/06 Chemicals");
        addParm("Dt", "Dt", "0.01");
        addParm("Converge Threshold", "Converge Threshold", "0.1");
        addParm("Verbose", "Verbose", "True",
                QStringList() << "True"
                              << "False");
        addParm("Debug Steps", "Debug Steps", "10");

        addParm("Auxin", "Auxin", "");
        addParm("Auxin Source", "Auxin Source", "0");
        addParm("Auxin Intercellular Diffusion", "Auxin Intercellular Diffusion", "1");
        addParm("Auxin Cell Permeability", "Auxin Cell Permeability", "0.2");
        addParm("Auxin Basal Production Rate", "Auxin Basal Production Rate", "0");
        addParm("Auxin QC Basal Production Rate", "Auxin QC Basal Production Rate", "0");
        addParm("Auxin SCN Basal Production Rate", "Auxin SCN Basal Production Rate", "0");
        addParm("Auxin Decay Rate", "Auxin Decay Rate", "0.0125");
        addParm("Auxin Max Decay Rate", "Auxin Decay Rate", "0.125");
        addParm("Auxin Max Amount Cell", "Auxin Max Amount (auxin per nm squared)", "3");
        addParm("Auxin Max Amount Edge", "Auxin Max Amount (auxin per nm)", "10");
        addParm("Pin1", "Pin1", "");
        addParm("Pin1 Basal Production Rate", "Pin1 Basal Production Rate", "0.2");
        addParm("Pin1 Max Auxin-induced Expression", "Pin1 Max Auxin-induced Expression", "30");
        addParm("Pin1 Half-max Auxin-induced K", "Pin1 Half-max Auxin-induced K", "0.05");
        addParm("Pin1 Max Concentration", "Pin1 Concentration", "2");
        //addParm("Pin1 Half-max Auxin-induced Decay", "Pin1 Half-max Auxin-induced Decay", "2000");
        //addParm("Pin1 Max Auxin-induced Decay", "Pin1 Max Auxin-induced Decay", "1");
        addParm("Pin1 Cytoplasmic Decay Rate", "Pin1 Cytoplasmic Decay Rate", "0.08");
        addParm("Pin1 Membrane Max Decay Rate", "Pin1 Membrane Max Decay Rate", "0.8");
        addParm("Pin1 Max Trafficking Rate", "Pin1 Max Trafficking Rate", "1");
        addParm("Pin1 Max Amount Edge", "Pin1 Max Amount Edge( auxin per nm)", "15");
        addParm("Pin1-auxin export rate", "Pin1-auxin export rate", "1.4");
        addParm("Pin1 Sensitivity Suppression by Auxin Amount",
                "Pin1 Sensitivity Suppression by Auxin Amount (auxin per nm squared)", "400"); //// be careful
        addParm("Pin1 Sensitivity Suppression by Auxin Max Cell",
                "Pin1 Sensitivity Suppression by Auxin Max Cell (auxin per nm squared)", "True", QStringList() << "True" << "False" ); //// be careful
        addParm("Simulate PIN4", "Simulate PIN4", "False", QStringList() << "True" << "False" );
        addParm("Columella Auto-Efflux", "Columella Auto-Efflux", "True", QStringList() << "True" << "False" );
        addParm("Pin1 Sensitivity MF K", "Pin1 Sensitivity MF K", "0");
        addParm("Pin1 Sensitivity Auxin-flux K", "Pin1 Sensitivity Auxin-flux K", "3");
        addParm("Pin1 Sensitivity Geometry K", "Pin1 Sensitivity Geometry K", "2");
        addParm("Pin1 Sensitivity MF+Auxin-flux K", "Pin1 Sensitivity MF+Auxin-flux K", "3");
        addParm("Pin1 Sensitivity MF+Geometry K", "Pin1 Sensitivity MF+Geometry K", "0");
        addParm("Pin1 Sensitivity Auxin-flux+Geometry K", "Pin1 Sensitivity Auxin-flux+Geometry K", "0");
        addParm("Pin1 Sensitivity Average Method", "Pin1 Sensitivity Average Method", "Soft-max",
                                                        QStringList() << "Soft-max"
                                                                      << "Arithmetic Average");

        addParm("Auxin Polarity Method", "Auxin Polarity Method", "Flow",
                QStringList() << "PINOID"
                              << "Flow");
        addParm("Auxin-Flux Impact Half-max", "Auxin-Flux Impact Half-max", "0.1");
        addParm("PINOID Impact Half-max", "PINOID Impact Half-max", "0.1");
        addParm("MF Impact Half-max",
                "MF Impact Half-max", "0.5");
        addParm("Geometry Impact Half-max",
                "Geometry Impact Half-max", "0.5");
        addParm("Aux1", "Aux1", "");
        addParm("Aux1 Basal Production Rate", "Aux1 Basal Production Rate", "1");
        addParm("Aux1 Half-max Auxin-induced K", "Aux1 Half-max Auxin-induced K", "0.01");
        addParm("Aux1 Max Auxin-induced Expression", "Aux1 Max Auxin-induced Expression", "30");
        addParm("Aux1 Max Concentration", "Aux1 Max Concentration", "2");
        addParm("Aux1 Cytoplasmic Decay Rate", "Aux1 Cytoplasmic Decay Rate", "0.08");
        addParm("Aux1 Max Trafficking Rate", "Aux1 Max Trafficking Rate", "1");
        addParm("Aux1 Max Amount Edge", "Aux1 Max Amount Edge", "15");
        addParm("Aux1-auxin import rate", "Aux1-auxin import rate", "1");
        addParm("Division Inhibitor", "Division Inhibitor", "");
        addParm("Division Inhibitor Basal Production Rate", "Division Inhibitor Basal Production Rate", "0");
        addParm("Division Inhibitor Max Promoter-induced Expression", "Division Inhibitor Max Promoter-induced Expression", "20");
        addParm("Division Inhibitor Half-max Promoter-induced K", "Division Inhibitor Half-max Promoter-induced K", "5"); // 20 or more for the mutant
        addParm("Division Inhibitor Half-max Promoter-induced n", "Division Inhibitor Half-max Promoter-induced n", "2");
        addParm("Division Inhibitor Decay Rate", "Division Inhibitor Decay Rate", "0.01");
        addParm("Division Inhibitor Permeability", "Division Inhibitor Permeability", "0");
        addParm("Division Promoter", "Division Promoter", "");
        addParm("Division Promoter Basal Production Rate", "Division Promoter Basal Production Rate", "0");
        addParm("Division Promoter Max Auxin-induced Expression", "Division Promoter Max Auxin-induced Expression", "20");
        addParm("Division Promoter Half-max Auxin-induced K", "Division Promoter Half-max Auxin-induced K", "2");
        addParm("Division Promoter Half-max Auxin-induced n", "Division Promoter Half-max Auxin-induced n", "4");
        addParm("Division Promoter Decay Rate", "Division Promoter Decay Rate", "0.01");
        addParm("Division Promoter Permeability", "Division Promoter Permeability", "1"); // 1 for the data, 5 for the figure
        addParm("Phosphorilation", "Phosphorilation", "");
        addParm("PINOID Basal Production Rate", "PINOID Basal Production Rate", "10");
        addParm("PP2A Basal Production Rate", "PP2A Basal Production Rate", "10");
        addParm("PINOID Dilution Rate", "PINOID Dilution Rate", "0.1");
        addParm("PP2A Dilution Rate", "PP2A Dilution Rate", "1");
        addParm("PINOID Decay Rate", "PINOID Decay Rate", "0.08");
        addParm("PP2A Decay Rate", "PP2A Decay Rate", "0.08");
        addParm("PINOID Trafficking Rate", "PINOID Trafficking Rate", "1");
        addParm("PP2A Trafficking Rate", "PP2A Trafficking Rate", "0.01");
        addParm("PINOID Max Amount Edge", "PINOID Max Amount Edge", "20");
        addParm("PP2A Max Amount Edge", "PP2A Max Amount Edge", "20");
        addParm("PINOID Displacing K", "PINOID Displacing K", "10");
        addParm("PINOID Fluidity K", "PINOID Fluidity K", "0.1");
        addParm("Auxin-PP2A Relief T", "Auxin-PP2A Relief T", "1");
        addParm("Auxin-PP2A Relief K", "Auxin-PP2A Relief K", "1");
        addParm("PINOID-PP2A Trafficking Toggle K", "PINOID-PP2A Trafficking Toggle K", "100");
        addParm("PINOID-PP2A Disassociation Toggle K", "PINOID-PP2A Disassociation Toggle K", "100");
        addParm("Geom-PP2A Relief T", "Geom-PP2A Relief T", "0.1");
        addParm("Geom-PP2A Relief K", "Geom-PP2A Relief K", "0.5");
        addParm("MF-PP2A Relief T", "Geom-PP2A Relief T", "0.1");
        addParm("MF-PP2A Relief K", "Geom-PP2A Relief K", "0.5");
        addParm("Tissue Process", "Name of Tissue Process", "Model/Root/03 Cell Tissue");
        addParm("Chemicals Solver Process",
                "Name of Chemicals solver process",
                "Model/Root/061 Chemicals Solver");
    }
    bool initialize(QWidget* parent);
    bool rewind(QWidget* parent);
    bool step();
    bool update();
    double calcNorm();
    Point8d calcDerivatives(CCIndex v, const CCStructure& csDual, const CCIndexDataAttr& indexAttr, Tissue::CellDataAttr &cellAttr,
                            Tissue::FaceDataAttr &faceAttr,
                            Tissue::EdgeDataAttr& edgeAttr, Tissue::VertexDataAttr& vMAttr);
    void calcDerivsEdge(const CCStructure &cs, const CCIndexDataAttr &indexAttr,
                        Tissue::EdgeDataAttr& edgeAttr,
                        CCIndex e, double Dt);
    void calcDerivsCell(const CCStructure& cs,
                           const CCStructure& csDual,
                           const CCIndexDataAttr& indexAttr,
                           Tissue::CellDataAttr& cellAttr,
                           Tissue::EdgeDataAttr& edgeAttr,
                           int label, double Dt);
    tbb::interface5::concurrent_unordered_map<string, double> debugs;
private:
    double Dt = 0;
    Tissue* tissueProcess = 0;
    Mesh* mesh = 0;
    CCIndexDataAttr* indexAttr = 0;
    Tissue::CellDataAttr* cellAttr = 0;
    Tissue::FaceDataAttr* faceAttr = 0;
    Tissue::VertexDataAttr* vMAttr = 0;
    double  userTime = 0;
    double convergeThresh = 0.01;
    int debug_step = 0;
};


struct MDXSubdivideX : public Subdivide
{
  MDXSubdivideX() {}
  MDXSubdivideX(Mesh &mesh);

  void splitCellUpdate(Dimension dim, const CCStructure &cs, const CCStructure::SplitStruct &ss,
      CCIndex otherP = CCIndex(), CCIndex otherN = CCIndex(), double interpPos = 0.5);

  Mesh *mesh = 0;
  CCIndexDataAttr *indexAttr = 0;
  Subdivide *nextSubdivide = 0;

  std::vector<CCIndexDoubleAttr*> attrVec;
};


class CellDivision : public Process, virtual public Cell2dDivideParms {
public:
    CellDivision(const Process& process)
        : Process(process) {
        //setName("03 Divide Cells");
        setDesc("Divide cells Parameters for cell division on tissue mesh.");
        setIcon(QIcon(":/images/CellDivide.png"));
        addParm("Verbose", "Verbose", "True",
                QStringList() << "True"
                              << "False");
        addParm("Division Algorithm",
                "Name of the Division Algorithm",
                "1",
                QStringList() << "0"
                              << "1"
                              << "2"  );
        addParm("Max Joining Distance",
                "When using division algorithm 2, the closest existing division point will be used as \
                joining point, up to this maximum distance",
                "1" );
        addParm("Minimum Polarity Vector Norm",
                "Minimum Polarity Vector Norm",
                "0.05" );
        addParm("Max Division Time",
                "Max Division Time",
                "50");
        addParm("Min Division Time",
                "Min Division Time",
                "2");
        addParm("Division Meristem Size",
                "Division Meristem Size",
                "200");
        addParm("Division Promoter Level",
                "Division Promoter Level",
                "0.05" );
        addParm("Division half-probability by Cell Size Ratio",
                "Division half-probability by Cell Size Ratio",
                "1.0" );
        addParm("Division half-probability by Inhibitor",
                "Division half-probability by Inhibitor",
                "0.03" );
        addParm("Division Control",
                "Division Control",
                "False",
                                QStringList() << "True"
                                              << "False");
        addParm("Ignore Cell Type",
                "Ignore Cell Type",
                "False",
                                QStringList() << "True"
                                              << "False");
        addParm("Root Process", "Name of the process for the Root", "Model/Root/01 Root");
        addParm("Remesh", "Remesh", "Model/Root/02 Remesh");
        addParm("Triangulate Faces Process", "Triangulate Faces", "Model/Root/Triangulate Faces");
        addParm("ClearCells Process", "ClearCells", "Model/Root/ClearCells");
        addParm("Tissue Process", "Name of Tissue Process", "Model/Root/03 Cell Tissue");

    }

    ~CellDivision() {}

    // Initialize simulation, called from GUI thread
    bool initialize(QWidget* parent);

    // Process the parameters
    bool processParms();

    // Run a step of cell division
    bool step() {
        Mesh* mesh = currentMesh();
        return step(mesh, &subdiv);
    }

    // Run a step of cell division
    virtual bool step(Mesh* mesh, Subdivide* subdiv);

    // Subdivide object
    MDXSubdivideX* subdivider() {
        return &subdiv;
    }
    Tissue::FaceDataAttr* FaceDataAttr = 0;

private:
    MDXSubdivideX subdiv;
    bool Verbose = false;
    Root* rootProcess = 0;
    Remesh* remeshProcess = 0;
    ClearCells* clearCellsProcess = 0;
    TriangulateFacesX* triangulateProcess = 0;
    Tissue* tissueProcess = 0;
    Point3d divVector;
};




class Splitter : virtual public mdx::Subdivide {
public:
    Splitter(bool cellDivision = false) {
        this->cellDivision = cellDivision;
    }

    void splitCellUpdate(Dimension dim,
                         const CCStructure& cs,
                         const CCStructure::SplitStruct& ss,
                         CCIndex otherP = CCIndex(),
                         CCIndex otherN = CCIndex(),
                         double interpPos = 0.5);

    MDXSubdivideX mdx;
    Tissue::Subdivide mechanics;
    bool cellDivision = false;
};

class RootDivide : public CellDivision {
public:
    RootDivide(const Process& process)
        : CellDivision(process) {
        setName("Model/Root/04 Divide Cells");

        addParm("Cell Division Enabled", "Cell Division Enabled", "True",
                                                        QStringList() << "False"
                                                                      << "True");
        addParm("Manual Cell Division Enabled", "Manual Cell Division Enabled", "False",
                                                        QStringList() << "False"
                                                                      << "True");
        setParmDefault("Cell Max Area", "100");
        setParmDefault("Cell Wall Min", "0.1");
    }

    // Initialize to grab subdivider
    bool initialize(QWidget* parent);

    // Run a step of cell division
    bool step(double Dt);

    Splitter& subdivider() {
        return subdiv;
    }

private:
    Splitter subdiv;

};


/**
 * \class TriangulateFaces CellMakerProcesses2D.hpp <CellMakerProcesses2D.hpp>
 * takes all faces of the CC and triangulates them using triangle
 *
 * TriangulateFaces
 */
class TriangulateFacesX : public Process {
public:
    TriangulateFacesX(const Process& process)
        : Process(process) {
        setName("Model/Root/Triangulate Faces");
        setDesc("Triangulate all polygonal faces.");
        addParm("Max Area", "Max area for triangles", "100");
        addParm("Destructive", "Destructive", "False",
                                                        QStringList() << "False"
                                                                      << "True");

    }

    bool triangulateFace(CCIndex f, CCStructure& cs,
                         CCIndexDataAttr& indexAttr,
                         double maxArea);

    bool step() {
        Mesh* mesh = getMesh("Mesh 1");
        if(!mesh or mesh->file().isEmpty())
            throw(QString("Root::TriangulateFacesX No current mesh"));

        CCIndexDataAttr& indexAttr = mesh->indexAttr();
        CCStructure& cs = mesh->ccStructure("Tissue");
        for(CCIndex f : cs.faces())
            if(indexAttr[f].selected) {
                triangulateFace(f, cs, indexAttr, parm("Max Area").toDouble());
                step(false);
                return false;
            }
        return step(parm("Destructive") == "True");
    }


    bool step(bool destructive) {
        Mesh* mesh = currentMesh();
        if(!mesh)
            throw QString("%1 No current mesh").arg(name());

        QString ccName = mesh->ccName();
        CCStructure& cs = mesh->ccStructure(ccName);
        CCIndexDataAttr& indexAttr = mesh->indexAttr();

        QString ccNameOut = ccName;
        CCStructure& csOut = mesh->ccStructure(ccNameOut);
        bool vVisible = mesh->drawParms(ccName).isGroupVisible("Vertices");
        bool eVisible = mesh->drawParms(ccName).isGroupVisible("Edges");
        bool fVisible = mesh->drawParms(ccName).isGroupVisible("Faces");
        bool result;
        if(destructive)
            result = step_destructive(cs, csOut, indexAttr, parm("Max Area").toDouble());
        else
            result = step_nondestructive(cs, indexAttr, parm("Max Area").toDouble());
        if(result)
            ccUpdateDrawParms(*mesh, ccName, ccNameOut);
        else {
            mesh->erase(ccNameOut);
            throw QString("%1 Triangulation of mesh failed").arg(name());
        }
        mesh->drawParms(ccNameOut).setGroupVisible("Vertices", vVisible);
        mesh->drawParms(ccNameOut).setGroupVisible("Edges", eVisible);
        mesh->drawParms(ccNameOut).setGroupVisible("Faces", fVisible);
        mesh->updateAll(ccNameOut);
        return false;
    }
    bool step_destructive(CCStructure& cs, CCStructure& csOut, CCIndexDataAttr& indexAttr, double maxArea);
    bool step_nondestructive(CCStructure& cs, CCIndexDataAttr& indexAttr, double maxArea);
};


class ClearCells : public Process {
public:
    ClearCells(const Process& process)
        : Process(process) {
        setName("Model/Root/ClearCells");
        setDesc("ClearCells.");
        addParm("AddFace", "AddFace", "Model/Root/Add Face");
    }


    CCIndex clearCell(int label);
    bool clearCell_old(int label);



    bool clearAllCells();

    bool step() {
        Mesh* mesh = getMesh("Mesh 1");
        if(!mesh or mesh->file().isEmpty())
            throw(QString("Root::ClearCells No current mesh"));
        CCIndexDataAttr& indexAttr = mesh->indexAttr();
        CCStructure& cs = mesh->ccStructure("Tissue");
        for(CCIndex f : cs.faces())
            if(indexAttr[f].selected) {
                clearCell(indexAttr[f].label);
                return false;
            }
        return clearAllCells();
    }

    AddFace *addFaceProcess = 0;

};

// This process is a big mess, it should be rewritten completely
// remesh the CS, a very messy job, remember to reinitialize
// all the processes as the CS is destroyed
class Remesh : public Process {
public:
    Remesh(const Process& process)
        : Process(process) {
        setName("Model/Root/02 Remesh");
        setDesc("Remesh");
        setIcon(QIcon(":/images/Recycle.png"));
        addParm("Destructive", "Destructive", "False",
                                                        QStringList() << "False"
                << "True");
        addParm("Type", "Type", "Hard",
                                                        QStringList() << "Hard"
                << "Soft");
        addParm("Split Edges Max Length", "Split Edges Max Length", "3");

        addParm("Remeshing Min Area", "Minimum triangle area that triggers remeshing", "0");
        addParm("Remeshing Max Area", "Maximum triangle area that triggers remeshing", "100");
        addParm("Tissue Process", "Name of Tissue Process", "Model/Root/03 Cell Tissue");
        addParm("Triangulate Faces", "Triangulate Faces", "Model/Root/Triangulate Faces");
        addParm("ClearCells", "ClearCells", "Model/Root/ClearCells");
        addParm("SplitEdges", "SplitEdges", "Mesh/Structure/Split Edges");
        addParm("Polygonize", "Polygonize", "Tools/Cell Maker/Mesh 2D/Tools/Polygonize Triangles");
    }

    bool initialize(QWidget* parent) {
        if(!getProcess(parm("Tissue Process"), tissueProcess))
            throw(QString("Root::initialize Cannot make Tissue Process:") + parm("TissueProcess"));
        if(!getProcess(parm("Triangulate Faces"), triangulateProcess))
            throw(QString("Root::initialize Cannot make Triangulate Faces") +
                  parm("Triangulate Faces"));
        if(!getProcess(parm("ClearCells"), clearCellsProcess))
            throw(QString("Root::initialize Cannot make ClearCells") + parm("ClearCells"));
        if(!getProcess(parm("SplitEdges"), splitEdgesProcess))
            throw(QString("Root::initialize Cannot make SplitEdges") + parm("SplitEdges"));
        if(!getProcess(parm("Polygonize"), polygonizeProcess))
            throw(QString("Root::initialize Cannot make Polygonize") + parm("Polygonize"));

        Mesh* mesh = getMesh("Mesh 1");
        indexAttr = &mesh->indexAttr();
        widget_parent = parent;
        return true;
    }

    void polygonize(bool destructive = false) {
        Mesh* mesh = getMesh("Mesh 1");
        CCStructure& cs = mesh->ccStructure("Tissue");
        if(destructive)
            polygonizeTriangles(cs, cs, *indexAttr, cs.faces(), 0);
        else
            clearCellsProcess->step();
    }


    void check() {
        Mesh* mesh = getMesh("Mesh 1");
        Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
        double maxArea = parm("Remeshing Max Area").toDouble();
        double minArea = parm("Remeshing Min Area").toDouble();
        std::set<int> cellToRemesh;

        // check if we should remesh
        for(auto c : cellAttr) {
            Tissue::CellData& cD = cellAttr[c.first];
            if(c.second.label == -1) // ?
                cellToRemesh.insert(cD.label);
            for(CCIndex f : *cD.cellFaces)
                if((*indexAttr)[f].measure > maxArea ||
                   (*indexAttr)[f].measure < minArea) {
                    mdxInfo << "Face " << f << " : "
                            << " label " << (*indexAttr)[f].label << " of size " << (*indexAttr)[f].measure
                            << " triggered remeshing" << endl;
                    cellToRemesh.insert(cD.label);
                }
        }

        if(cellToRemesh.size() == 0)
            return;
/*
        for(int label : cellToRemesh) {
            CCIndex ft = clearCellsProcess->clearCell(label);
            updateGeometry(cs, (*indexAttr));
            (*indexAttr)[ft].selected = true;
        }

        triangulateProcess->step();
        tissueProcess->initialize();*/
        step(true, false, false);

        return ;
    }

    bool step(bool forceRemesh=false, bool destructive = false, bool soft = false) {
        Mesh* mesh = getMesh("Mesh 1");
        CCStructure& cs = mesh->ccStructure("Tissue");

        mdxInfo << "Remeshing!" << (soft ? "(Soft)" : "(Hard)") << endl;
        Splitter subdiv;
        // Setup subdivision objects
        subdiv.mdx = MDXSubdivideX(*mesh);
        subdiv.mechanics =
            Tissue::Subdivide(*indexAttr,
                              mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData"),
                              mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData"),
                              mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData"),
                              mesh->attributes().attrMap<int, Tissue::CellData>("CellData"));

        // take a copy of current cs
        CCStructure csCurr = cs;
        // Set Z zero
        for(CCIndex v : cs.vertices())
            (*indexAttr)[v].pos[2] = 0;
        // create the tissue, so that the dual one is based on one-faced cells
        tissueProcess->initialize(widget_parent);
        // let's wipe delete cells' internal faces
        if(!soft) {
            mdxInfo << "Clearing Cells..." << endl;
            polygonize(destructive);
            splitEdgesProcess->run(*mesh, cs, parm("Split Edges Max Length").toDouble(), &subdiv);
            // create the tissue, so that the dual one is based on one-faced cells
            tissueProcess->initialize(widget_parent);
            // triangulate, this will erase the current cs and create a new one!
            mdxInfo << "Triangulating..." << endl;
            triangulateProcess->step(destructive);
        }
        // restore all tissue's properties, based on the old cs
        mdxInfo << "Restoring Cell Tissue..." << endl;
        tissueProcess->restore(csCurr);
        // update geometries as usual
        tissueProcess->step(EPS);
        // initialize the provided processes
        for(Process* p : processes)
            p->initialize(widget_parent);
        // check the CS for consistency
        if(!verifyCCStructure(cs, *indexAttr))
            throw(QString("Remesh::step The Cell Complex is broken"));

        return false;
    }

    bool run(){
        bool soft = parm("Type") == "Soft";
        bool destructive = parm("Destructive") == "True";
        step(true, destructive, soft);
        return false;
    }

    void setProcesses(std::vector<Process*> processes) {
        this->processes = processes;
    }

private:
    QWidget* widget_parent = 0;
    CCIndexDataAttr* indexAttr = 0;
    Tissue* tissueProcess = 0;
    ClearCells* clearCellsProcess = 0;
    mdx::SplitEdges* splitEdgesProcess = 0;
    Process* polygonizeProcess = 0;
    TriangulateFacesX* triangulateProcess = 0;
    std::vector<Process*> processes;
};


// Main model class
class Root : public Process {
public:
    Root(const Process& process)
        : Process(process) {
        setName("Model/Root/01 Root");
        addParm("Max Mechanical Iterations", "Max Mechanical Iterations", "1");
        addParm("Max Chemical Iterations", "Max Chemical Iterations", "1");
        addParm("Mechanics Enabled",
                "Mechanics Enabled",
                "True",
                QStringList() << "True"
                              << "False");
        addParm("Chemicals Enabled",
                "Chemicals Enabled",
                "True",
                QStringList() << "True"
                              << "False");
        addParm("Growth Enabled",
                "Growth Enabled",
                "True",
                QStringList() << "True"
                              << "False");
        addParm("Cell Division Enabled",
                "Cell Division Enabled",
                "True",
                QStringList() << "True"
                              << "False");
        addParm("Remesh at start",
                "Remesh at start",
                "Hard",
                QStringList() << "Hard"
                              << "Soft"
                              << "None");
        addParm("Remesh during execution",
                "Remesh during execution",
                "True",
                QStringList() << "True"
                              << "False");
        addParm("Mesh Update Timer", "Mesh Update Timer", "1");
        addParm("Snapshots Timer",
                "Time frames between snapshots",
                "0");
        addParm("Snapshots Directory",
                "Snapshots Directory",
                "");
        addParm("Debug",
                "Debug",
                "True",
                QStringList() << "False"
                              << "True");
        addParm("Debug File", "Debug File", "debug.csv");
        addParm("Frame fixed on QC", "Frame fixed on QC", "0");
        addParm("Frame fixed on Substrate", "Frame fixed on Substrate", "0");
        addParm("Execution Time", "Execution Time", "0");
        addParm("Output Mesh", "Output Mesh", "output.mdxm");
        addParm("Tissue Process", "Name of Tissue Process", "Model/Root/03 Cell Tissue");
        addParm("Mechanical Solver Process",
                "Name of Mechanical Solver Process",
                "Model/Root/051 Mechanical Solver");
        addParm("Mechanics Process", "Name of Mechanics Process", "Model/Root/05 Mechanics");
        addParm("Mechanical Growth Process",
                "Name of the process for Mechanical Growth",
                "Model/Root/052 Mechanical Growth");
        addParm(
            "Chemicals Process", "Name of the process for Chemicals", "Model/Root/06 Chemicals");
        addParm("Chemicals Solver Process",
                "Name of the process for Chemicals Solver",
                "Model/Root/061 Chemicals Solver");
        addParm("Divide Process",
                "Name of the process for Cell Division",
                "Model/Root/04 Divide Cells");
        addParm("Delete Cell Process",
                "Name of the process for Delete Cell",
                "Model/Root/31 Delete Cell");
        addParm("Execute Test Process",
                "Name of the process for Execute Test",
                "Model/Root/4 Execute Test");
        addParm("Set Global Attr Process",
                "Name of the process for Set Global Attr",
                "Model/Root/23 Set Global Attr");
        addParm("Triangulate Faces", "Triangulate Faces", "Model/Root/Triangulate Faces");
        addParm("Remesh", "Remesh", "Model/Root/02 Remesh");
        addParm("SplitEdges", "SplitEdges", "Mesh/Structure/Split Edges");
        addParm("SaveView", "SaveView", "Tools/System/Save View");
        addParm("SaveMesh", "SaveMesh", "Mesh/System/Save");
    }

    // Initialize the model
    bool initialize(QWidget* parent);

    // Run the model
    // bool run();
    bool step();

    bool finalize();

    // Rewind the model (reloads the mesh)
    bool rewind(QWidget* parent);

    Mesh* mesh = 0;
    QWidget* widget_parent = 0;
    MechanicalGrowth* mechanicalGrowthProcess = 0;
    Chemicals* chemicalsProcess = 0;
    Mechanics* mechanicsProcess = 0;
    Tissue* tissueProcess = 0;
    RootDivide* divideProcess = 0;
    DeleteCell* deleteProcess = 0;
    ExecuteTest* executeTestProcess = 0;
    SetGlobalAttr* setGlobalAttrProcess = 0;
    TriangulateFacesX* triangulateProcess = 0;
    Remesh* remeshProcess = 0;
    Process* splitEdgesProcess = 0;
    SaveViewFile* saveViewProcess = 0;
    MeshSave* saveMeshProcess = 0;
    std::vector<Process*> processes;
    double userTime = 0;

private:
    bool debugging = false;
    bool process_started = false;
    bool snapshot = false;
    string snapshotDir;
    int screenShotCount = 0;
    bool mechanicsEnabled = true;
    bool chemicalsEnabled = true;
    bool growthEnabled = true;
    bool divisionEnabled = true;
    int stepCount = 0, prevStepCount = 0;
    clock_t begin_clock, prev_clock;
    std::ofstream output_file;
    int maxMechanicsIter = 0, maxChemicalIter = 0;
    bool mechanicsProcessConverged = false;
};


class DeleteCell : public Process {
public:
    DeleteCell(const Process& process)
        : Process(process) {
        setName("Model/Root/31 Delete Cell");
        setDesc("DeleteCell.");
        addParm("Tissue Process", "Name of Tissue Process", "Model/Root/03 Cell Tissue");

    }

    bool step(std::set<int> labels) {
        Mesh* mesh = getMesh("Mesh 1");
        CCStructure& cs = mesh->ccStructure("Tissue");
        Tissue::CellDataAttr &cellAttr =
            mesh->attributes().attrMap<int, Tissue::CellData>(
                "CellData");

        if(!getProcess(parm("Tissue Process"), tissueProcess))
            throw(QString("Root::initialize Cannot make Tissue Process") + parm("Tissue Process"));

        std::set<int> to_delete;
        for(auto c : cellAttr) {
            Tissue::CellData& cD = cellAttr[c.first];
            if(labels.find(cD.label) != labels.end())
                 to_delete.insert(cD.label);
        }
        for(int label : to_delete)
                deleteCell(label);

        tissueProcess->initialize();
        tissueProcess->restore(cs);
        tissueProcess->step(EPS);
        mesh->updateAll();

        return false;
    }

    bool step() {
        Mesh* mesh = getMesh("Mesh 1");
        CCStructure& cs = mesh->ccStructure("Tissue");
        CCIndexDataAttr& indexAttr = mesh->indexAttr();
        if(!getProcess(parm("Tissue Process"), tissueProcess))
            throw(QString("Root::initialize Cannot make Tissue Process") + parm("Tissue Process"));

        std::set<int> to_delete;
        for(CCIndex f : cs.faces())
            if(indexAttr[f].selected)
                to_delete.insert(indexAttr[f].label);
        for(int label : to_delete)
                deleteCell(label);

        tissueProcess->initialize();
        tissueProcess->restore(cs);
        tissueProcess->step(EPS);
        mesh->updateAll();

        return false;
    }


    void deleteCell(int label) {
        Mesh* mesh = getMesh("Mesh 1");
        CCStructure& cs = mesh->ccStructure("Tissue");
        CCStructure& csVisual = mesh->ccStructure("TissueVisual");

        CCIndexDataAttr& indexAttr = mesh->indexAttr();
        Tissue::CellDataAttr &cellAttr =
            mesh->attributes().attrMap<int, Tissue::CellData>(
                "CellData");
        Tissue::FaceDataAttr &faceAttr =
            mesh->attributes().attrMap<CCIndex, Tissue::FaceData>(
                "FaceData");
        Tissue::VertexDataAttr &vMAttr =
            mesh->attributes().attrMap<CCIndex, Tissue::VertexData>(
                "VertexData");
        if(cellAttr.find(label) == cellAttr.end())
            throw(QString("DeleteCell: cell labelled" + QString::number(label) + " not found"));
        Tissue::CellData& cD = cellAttr[label];
        std::set<CCIndex> to_delete_faces;
        std::set<CCIndex> to_delete_edges;
        std::set<CCIndex> to_delete_vertices;
        for(CCIndex f : *cD.cellFaces) {
            std::set<int> labels;
            Tissue::FaceData fD = faceAttr[f];
            for(CCIndex e : cs.incidentCells(f, 1)) {
                labels.clear();
                for(CCIndex fn : cs.incidentCells(e, 2))
                    labels.insert(indexAttr[fn].label);
                if(labels.size() == 1)
                    to_delete_edges.insert(e);
            }
            for(CCIndex v : cs.incidentCells(f, 0)) {
                labels.clear();
                for(CCIndex fn : cs.incidentCells(v, 2))
                    labels.insert(indexAttr[fn].label);
                if (labels.size() == 1)
                    to_delete_vertices.insert(v);
            }
            if(csVisual.hasCell(fD.G_v0)) {
                csVisual.deleteCell(fD.G_e1);
                csVisual.deleteCell(fD.G_e2);
                csVisual.deleteCell(fD.G_v0);
                csVisual.deleteCell(fD.G_v1);
                csVisual.deleteCell(fD.G_v2);
            }
            if(csVisual.hasCell(fD.a1_v1)) {
                csVisual.deleteCell(fD.a1_e);
                csVisual.deleteCell(fD.a2_e);
                csVisual.deleteCell(fD.a1_v1);
                csVisual.deleteCell(fD.a1_v2);
                csVisual.deleteCell(fD.a2_v1);
                csVisual.deleteCell(fD.a2_v2);
            }
            to_delete_faces.insert(f);
        }
        for(auto i : to_delete_faces)
            cs.deleteCell(i);
        for(auto i : to_delete_edges)
            cs.deleteCell(i);
        for(auto v : to_delete_vertices) {
            Tissue::VertexData vD = vMAttr[v];
            if(csVisual.hasCell(vD.V_e)) {
                csVisual.deleteCell(vD.V_e);
                csVisual.deleteCell(vD.V_v0);
                csVisual.deleteCell(vD.V_v1);
            }
            cs.deleteCell(v);
        }
        if(csVisual.hasCell(cD.a1_e)) {
            csVisual.deleteCell(cD.a1_e);
            csVisual.deleteCell(cD.a1_v1);
            csVisual.deleteCell(cD.a1_v2);
            csVisual.deleteCell(cD.a2_e);
            csVisual.deleteCell(cD.a2_v1);
            csVisual.deleteCell(cD.a2_v2);
        }
        if(csVisual.hasCell(cD.auxinFlux_f)) {
            std::set<CCIndex> to_delete;
            for(CCIndex e : csVisual.incidentCells(cD.auxinFlux_f, 1))
                to_delete.insert(e);
            csVisual.deleteCell(cD.auxinFlux_f);
            for(CCIndex e : to_delete)
                csVisual.deleteCell(e);
            csVisual.deleteCell(cD.auxinFlux_e);
            csVisual.deleteCell(cD.auxinFlux_v1);
            csVisual.deleteCell(cD.auxinFlux_v2);
            csVisual.deleteCell(cD.auxinFlux_v3);
            csVisual.deleteCell(cD.auxinFlux_v4);
        }
        mdxInfo << "Cell " << cD.label << " in position " << label <<  " removed" << endl;
        cellAttr.erase(label);
    }

private:
    Tissue* tissueProcess = 0;

};



class ExecuteTest : public Process {
private:
    Tissue* tissueProcess = 0;
    DeleteCell* deleteProcess = 0;
    MechanicalGrowth * mechanicalGrowthProcess = 0;
    DeleteSelection* deleteSelectionProcess = 0;
    Remesh* remeshProcess = 0;
    std::map<int, double> cellsProdRates;
    double auxinConc = 0;
    int source_removal_count = 0, alternate_source_count = 0, auxin_overflow_count = 0, QC_ablation_count = 0, LRC_removal_count = 0,
        tip_removal_count = 0, pin2_removal_count = 0, pin1_removal_count = 0, aux1_removal_count = 0,aux1_overexpr_count = 0, strain_removal_count = 0, endodermal_cells_count = 0;
    bool PIN1_knockdown = false;
    bool AUX1_knockdown = false;
    bool AUX1_overexpr = false;
    bool LRC_removed = false;
    bool tip_removed = false;
    bool overflow = false;

public:
    ExecuteTest(const Process& process)
        : Process(process) {
        setName("Model/Root/4 Execute Test");
        setDesc("ExecuteTest.");
        addParm("QC Ablation", "QC Ablation", "0");
        addParm("Source Removal", "Source Removal", "0");
        addParm("Alternate Source Removal", "Alternate Source Removal", "0");
        addParm("Auxin Overflow", "Auxin Overflow", "0");
        addParm("TIP Removal Time", "TIP Removal Time", "0");
        addParm("LRC Removal Time", "LRC Removal Time", "0");
        addParm("PIN2 Knockdown Time", "PIN2 Knockdown Time", "0");
        addParm("PIN1 Knockdown Time", "PIN1 Knockdown Time", "0");
        addParm("AUX1 Knockdown Time", "AUX1 Knockdown Time", "0");
        addParm("AUX1 Overexpression Time", "AUX1 Overexpression Time", "0");
        addParm("Strain Removal Time", "Strain Removal Time", "0");
        addParm("Endodermal Cells Delete", "Endodermal Cells Delete", "0");
        addParm("Tissue Process", "Name of Tissue Process", "Model/Root/03 Cell Tissue");
        addParm("Delete Cell Process", "Delete Cell Process", "Model/Root/31 Delete Cell");
        addParm("Mechanical Growth Process",
                "Name of the process for Mechanical Growth",
                "Model/Root/052 Mechanical Growth");
        addParm("Remesh Process", "Remesh Process", "Model/Root/02 Remesh");
        addParm("Delete Selection Process", "Delete Selection Process", "Mesh/System/Delete Selection");

        source_removal_count = 0, alternate_source_count = 0, auxin_overflow_count = 0, QC_ablation_count = 0, LRC_removal_count = 0,
        tip_removal_count = 0, pin1_removal_count = 0, pin2_removal_count = 0, aux1_removal_count = 0,aux1_overexpr_count = 0, strain_removal_count = 0, endodermal_cells_count = 0;;
        PIN1_knockdown = false;
        AUX1_knockdown = false;
        AUX1_overexpr = false;
        LRC_removed = false;
        tip_removed = false;
        overflow = false;
        cellsProdRates.clear();
        auxinConc = 0;
    }


    bool rewind(QWidget* parent) {
        source_removal_count = 0, alternate_source_count = 0, auxin_overflow_count = 0, QC_ablation_count = 0, LRC_removal_count = 0,
        tip_removal_count = 0, pin2_removal_count = 0, pin1_removal_count = 0, aux1_removal_count = 0,aux1_overexpr_count = 0, strain_removal_count = 0, endodermal_cells_count = 0;
        PIN1_knockdown = false;
        AUX1_knockdown = false;
        AUX1_overexpr = false;
        LRC_removed = false;
        tip_removed = false;
        overflow = false;
        cellsProdRates.clear();
        auxinConc = 0;
        return false;
    }

    bool step() {return false;}

    bool step(int stepCount=0) {
        Mesh* mesh = getMesh("Mesh 1");
        CCStructure& cs = mesh->ccStructure("Tissue");
        CCIndexDataAttr& indexAttr = mesh->indexAttr();
        Tissue::EdgeDataAttr& edgeAttr =
            mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
        Tissue::CellDataAttr &cellAttr =
            mesh->attributes().attrMap<int, Tissue::CellData>(
                "CellData");
        Tissue::VertexDataAttr& vMAttr =
            mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");

        if(!getProcess(parm("Tissue Process"), tissueProcess))
            throw(QString("Root::initialize Cannot make Tissue Process") + parm("Tissue Process"));
        if(!getProcess(parm("Delete Cell Process"), deleteProcess))
            throw(QString("Root::initialize Cannot make Delete Cell Process:") + parm("Delete Cell Process"));
        if(!getProcess(parm("Mechanical Growth Process"), mechanicalGrowthProcess))
            throw(QString("Root::initialize Cannot make Mechanical Growth Process:") +
                  parm("Mechanical Growth Process"));
        if(!getProcess(parm("Remesh Process"), remeshProcess))
            throw(QString("Root::initialize Cannot make Remesh Process:") +
                  parm("Remesh Process"));
        if(!getProcess(parm("Delete Selection Process"), deleteSelectionProcess))
            throw(QString("Root::initialize Cannot make Delete Selection Process:") + parm("Delete Selection Process"));


        int QC_ablation_time = parm("QC Ablation").toInt();
        QC_ablation_count ++;
        if(QC_ablation_time > 0 && QC_ablation_count > QC_ablation_time) {
            QC_ablation_count = 0;
            for(auto c : cellAttr) {
                Tissue::CellData& cD = cellAttr[c.first];
                if(cD.type == Tissue::QC)
                     cD.type = Tissue::Undefined;
            }
        }

        int source_removal_time = parm("Source Removal").toDouble();
        source_removal_count ++;
        if(source_removal_time > 0 && source_removal_count > source_removal_time) {
            for(auto c : cellAttr) {
                Tissue::CellData& cD = cellAttr[c.first];
                if(cD.type == Tissue::Source)
                    cD.setType(Tissue::Substrate, vMAttr);
            }
        }

        int alternate_source_time = parm("Alternate Source Removal").toDouble();
        alternate_source_count ++;
        if(alternate_source_time > 0 && alternate_source_count > alternate_source_time) {
            alternate_source_count = 0;
            for(auto c : cellAttr) {
                Tissue::CellData& cD = cellAttr[c.first];
                if(cD.auxinProdRate > 0) {
                    cellsProdRates[cD.label] = cD.auxinProdRate;
                    cD.auxinProdRate = 0;
                }
                else if (cD.auxinProdRate == 0) {
                    cD.auxinProdRate = cellsProdRates[cD.label];
                }
            }
        }

        QStringList list_overflow = parm("Auxin Overflow").split(QRegExp(","));
        if(list_overflow.size() != 4)
             throw(QString("Specify the Auxin Overflow test as start,inverval,amout"));
        int start = list_overflow[0].toInt();
        int interval1 = list_overflow[1].toInt();
        int interval2 = list_overflow[2].toInt();
        int amount = list_overflow[3].toInt();
        if(stepCount > start) {
            auxin_overflow_count ++;
            if(!overflow && interval1 > 0 && auxin_overflow_count > interval1) {
                auxinConc = 0;
                for(auto c : cellAttr) {
                    Tissue::CellData& cD = cellAttr[c.first];
                    if(cD.type != Tissue::Source) {
                        cellsProdRates[cD.label] = cD.auxinProdRate;
                        cD.auxinProdRate = amount;
                        auxinConc += cD.auxin;
                    }
                }
                overflow = true;
                auxin_overflow_count = 0;
                auxinConc /= cellAttr.size();

            }
            if(overflow && interval2 > 0 && auxin_overflow_count > interval2) {
                for(auto c : cellAttr) {
                    Tissue::CellData& cD = cellAttr[c.first];
                    if (cD.type != Tissue::Source) {
                        cD.auxinProdRate = cellsProdRates[cD.label];
                        cD.auxin = auxinConc;
                    }
                }
                overflow = false;
                auxin_overflow_count = 0;

            }

        }

        int tip_removal_time = parm("TIP Removal Time").toDouble();
        tip_removal_count ++;
        if(!tip_removed && tip_removal_time > 0 && tip_removal_count > tip_removal_time) {
            std::set<int> labels;
            Point3d VIcm;
            int VIcount = 0;
            for(auto c : cellAttr) {
                Tissue::CellData& cD = cellAttr[c.first];
                if(cD.type == Tissue::VascularInitial) {
                    VIcm += cD.centroid;
                    VIcount++;
                }
            }
            VIcm /= VIcount;
            /*
            for(auto c : cellAttr) {
                Tissue::CellData& cD = cellAttr[c.first];
                if(cD.type == Tissue::QC || cD.centroid.y() < QCcm.y())
                    labels.insert(cD.label);
            }
            deleteProcess->step(labels);
            */
            for(CCIndex v : cs.vertices())
                if(indexAttr[v].pos[1] < VIcm.y())
                    indexAttr[v].selected = true;
            deleteSelectionProcess->run(cs, indexAttr);
            remeshProcess->step(true, false, false);
            tip_removed = true;
        }

        int LRC_removal_time = parm("LRC Removal Time").toDouble();
        LRC_removal_count ++;
        if(LRC_removal_time > 0 && LRC_removal_count > LRC_removal_time && !LRC_removed) {
            std::set<int> labels;
            for(auto c : cellAttr) {
                Tissue::CellData& cD = cellAttr[c.first];
                if(cD.type == Tissue::LRC)
                    labels.insert(cD.label);
            }
            deleteProcess->step(labels);
            LRC_removed = true;
        }

        int pin2_removal_time = parm("PIN2 Knockdown Time").toDouble();
        pin2_removal_count ++;
        if(pin2_removal_time > 0 && pin2_removal_count > pin2_removal_time) {
            for(auto c : cellAttr) {
                Tissue::CellData& cD = cellAttr[c.first];
                if(cD.type == Tissue::LRC || cD.type == Tissue::Cortex || cD.type == Tissue::Epidermis || cD.type == Tissue::EpLrcInitial)
                    cD.Pin1 = 1;

            }
        }

        int pin1_removal_time = parm("PIN1 Knockdown Time").toDouble();
        pin1_removal_count ++;
        if(!PIN1_knockdown && pin1_removal_time > 0 && pin1_removal_count > pin1_removal_time) {
            for(auto c : cellAttr) {
                Tissue::CellData& cD = cellAttr[c.first];
                if(cD.type == Tissue::Vascular || cD.type == Tissue::VascularInitial || cD.type == Tissue::Pericycle || cD.type == Tissue::Endodermis) {
                    cD.pinProdRate /= 10;
                    cD.pinInducedRate /= 10;
                }
                PIN1_knockdown = true;
            }
        }

        int aux1_removal_time = parm("AUX1 Knockdown Time").toDouble();
        aux1_removal_count ++;
        if(!AUX1_knockdown && aux1_removal_time > 0 && aux1_removal_count > aux1_removal_time) {
            for(auto c : cellAttr) {
                Tissue::CellData& cD = cellAttr[c.first];
                if(cD.type != Tissue::Source && cD.type != Tissue::Substrate ) {
                    cD.Aux1 = 0;
                    cD.aux1InducedRate = 0;
                    //cD.aux1ProdRate /= 10;
                    for(CCIndex e : cD.perimeterEdges) {
                        Tissue::EdgeData& eD = edgeAttr[e];
                        eD.Aux1[cD.label] /= 10;
                    }
                }
                AUX1_knockdown = true;
           }

        }

        int aux1_overexpr_time = parm("AUX1 Overexpression Time").toDouble();
        aux1_overexpr_count ++;
        if(!AUX1_overexpr && aux1_overexpr_time > 0 && aux1_overexpr_count > aux1_overexpr_time) {
            for(auto c : cellAttr) {
                Tissue::CellData& cD = cellAttr[c.first];
                if(cD.type != Tissue::Source && cD.type != Tissue::Substrate ) {
                    cD.aux1MaxEdge *= 2;
                    cD.aux1ProdRate = 100;
                }
            }
            AUX1_overexpr = true;
        }


        int strain_removal_time = parm("Strain Removal Time").toDouble();
        strain_removal_count ++;
        if(strain_removal_time > 0 && strain_removal_count > strain_removal_time) {
               mechanicalGrowthProcess->setParm("MF Degradation", "0.5");
        }

        int endodermal_cells_time = parm("Endodermal Cells Delete").toDouble();
        endodermal_cells_count ++;
        if(stepCount > 1200 && endodermal_cells_time > 0 && endodermal_cells_count > endodermal_cells_time) {
            endodermal_cells_count = 0;
            std::vector<int> endos ;
            for(auto c : cellAttr) {
                Tissue::CellData& cD = cellAttr[c.first];
                if(cD.type == Tissue::Endodermis)
                    endos.push_back(cD.label);
            }
            // choose a random endodermal cell to delete
            std::set<int> to_delete;
            std::random_shuffle(endos.begin(), endos.end());
            to_delete.insert(*endos.begin());
            deleteProcess->step(to_delete);
        }

        return false;
    }


};



class PrintCellAttr : public Process {
public:
    PrintCellAttr(const Process& process)
        : Process(process) {
        setName("Model/Root/21 Print Cell Attr");
        setDesc("Print the cell type.");
    }
    bool step();
};


class DumpSignalInfo : public Process {
public:
    DumpSignalInfo(const Process& process)
        : Process(process) {
        setName("Model/Root/26 Dump Signal Info");
        setDesc("Dump Signal Info.");
        addParm("Signal",
                "Signal",
                "Chems: Auxin By Area");
    }

    bool step();
};





class SetCellAttr : public Process {
public:
    SetCellAttr(const Process& process)
        : Process(process) {
        setName("Model/Root/22 Set Cell Attr");
        addParm("Cell Type",
                "Cell Type",
                "QC",
                QStringList()
                              << "Undefined"
                              << "QC"
                              << "Columella"
                              << "ColumellaInitial"
                              << "CEI"
                              << "CEID"
                              << "Cortex"
                              << "Endodermis"
                            << "VascularInitial"
                            << "Vascular"
                              << "Pericycle"
                              << "EpLrcInitial"
                              << "Epidermis"
                              << "LRC"
                              << "Columella"
                              << "Substrate"
                              << "Source");
        addParm("Vertices Masses", "Vertices Masses", "1");
        addParm("Auxin production rate", "Auxin production rate", "0");
        addParm("Division Algorithm", "Division Algorithm", "-1");
        addParm("Periclinal Division", "Periclinal Division","",
                QStringList() << ""
                              << "False"
                              << "True");
        addParm("MF reorientation rate", "MF reorientation rate", "-1");
        addParm("Microfibril 1", "Microfibril 1", "0,1,0");
        addParm("Microfibril 2", "Microfibril 2", "0,-1,0");

    }
    bool step();
};


class SetGlobalAttr : public Process {
public:
    SetGlobalAttr(const Process& process)
        : Process(process) {
        setName("Model/Root/23 Set Global Attr");
        addParm("Mechanics Process", "Name of Mechanics Process", "Model/Root/05 Mechanics");
        addParm("Mechanical Growth Process",
                "Name of the process for Mechanical Growth",
                "Model/Root/052 Mechanical Growth");
        addParm("Divide Process", "Name of Divide Process", "Model/Root/04 Divide Cells");
        addParm("Undefined Wall EK", "", "-1");
        addParm("Undefined Wall CK", "", "-1");
        addParm("Undefined Shear EK", "", "-1");
        addParm("Undefined Shear CK", "", "-1");
        addParm("Undefined Turgor Pressure", "", "-1");
        addParm("Undefined Growth Factor", "", "0");
        addParm("Undefined Max area", "", "10000");
        addParm("Undefined MF reorientation rate", "", "0");
        addParm("QC Wall EK", "", "1");
        addParm("QC Wall CK", "", "1");
        addParm("QC Shear EK", "", "1");
        addParm("QC Shear CK", "", "1");
        addParm("QC Turgor Pressure", "", "0");
        addParm("QC Growth Factor", "", "0");
        addParm("QC Max area", "", "1000");
        addParm("QC MF reorientation rate", "", "0");
        addParm("ColumellaInitial Wall EK", "", "-1");
        addParm("ColumellaInitial Wall CK", "", "-1");
        addParm("ColumellaInitial Shear EK", "", "-1");
        addParm("ColumellaInitial Shear CK", "", "-1");
        addParm("ColumellaInitial Growth Factor", "", "1");
        addParm("ColumellaInitial Turgor Pressure", "", "1");
        addParm("ColumellaInitial Max area", "", "100");
        addParm("ColumellaInitial MF reorientation rate", "", "0");
        addParm("Columella Wall EK", "", "-1");
        addParm("Columella Wall CK", "", "-1");
        addParm("Columella Shear EK", "", "-1");
        addParm("Columella Shear CK", "", "-1");
        addParm("Columella Growth Factor", "", "1");
        addParm("Columella Turgor Pressure", "", "1");
        addParm("Columella MF reorientation rate", "", "0");
        addParm("Columella Max area", "", "150");
        addParm("VascularInitial Wall EK", "", "-1");
        addParm("VascularInitial Wall CK", "", "-1");
        addParm("VascularInitial Shear EK", "", "-1");
        addParm("VascularInitial Shear CK", "", "-1");
        addParm("VascularInitial Turgor Pressure", "", "-1");
        addParm("VascularInitial Growth Factor", "", "5");
        addParm("VascularInitial MF reorientation rate", "", "-1");
        addParm("VascularInitial Max area", "", "70");
        addParm("Vascular Wall EK", "", "-1");
        addParm("Vascular Wall CK", "", "-1");
        addParm("Vascular Shear EK", "", "-1");
        addParm("Vascular Shear CK", "", "-1");
        addParm("Vascular Turgor Pressure", "", "-1");
        addParm("Vascular Growth Factor", "", "5");
        addParm("Vascular MF reorientation rate", "", "0");
        addParm("Vascular Max area", "", "100");
        addParm("Pericycle Wall EK", "", "-1");
        addParm("Pericycle Wall CK", "", "-1");
        addParm("Pericycle Shear EK", "", "-1");
        addParm("Pericycle Shear CK", "", "-1");
        addParm("Pericycle Turgor Pressure", "", "-1");
        addParm("Pericycle Growth Factor", "", "1");
        addParm("Pericycle MF reorientation rate", "", "-1");
        addParm("Pericycle Max area", "", "75");
        addParm("Cortex Wall EK", "", "-1");
        addParm("Cortex Wall CK", "", "-1");
        addParm("Cortex Shear EK", "", "-1");
        addParm("Cortex Shear CK", "", "-1");
        addParm("Cortex Turgor Pressure", "", "-1");
        addParm("Cortex Growth Factor", "", "1");
        addParm("Cortex Max area", "", "75");
        addParm("Cortex MF reorientation rate", "", "-1");
        addParm("Endodermis Wall EK", "", "-1");
        addParm("Endodermis Wall CK", "", "-1");
        addParm("Endodermis Shear EK", "", "-1");
        addParm("Endodermis Shear CK", "", "-1");
        addParm("Endodermis Turgor Pressure", "", "-1");
        addParm("Endodermis Growth Factor", "", "1");
        addParm("Endodermis Max area", "", "50");
        addParm("Endodermis MF reorientation rate", "", "-1");
        addParm("Epidermis Wall EK", "", "-1");
        addParm("Epidermis Wall CK", "", "-1");
        addParm("Epidermis Shear EK", "", "-1");
        addParm("Epidermis Shear CK", "", "-1");
        addParm("Epidermis Turgor Pressure", "", "-1");
        addParm("Epidermis Growth Factor", "", "1");
        addParm("Epidermis Max area", "", "75");
        addParm("Epidermis MF reorientation rate", "", "-1");
        addParm("CEI Wall EK", "", "-1");
        addParm("CEI Wall CK", "", "-1");
        addParm("CEI Shear EK", "", "-1");
        addParm("CEI Shear CK", "", "-1");
        addParm("CEI Turgor Pressure", "", "-1");
        addParm("CEI Growth Factor", "", "1");
        addParm("CEI Max area", "", "150");
        addParm("CEI MF reorientation rate", "", "-1");
        addParm("CEID Wall EK", "", "-1");
        addParm("CEID Wall CK", "", "-1");
        addParm("CEID Shear EK", "", "-1");
        addParm("CEID Shear CK", "", "-1");
        addParm("CEID Turgor Pressure", "", "-1");
        addParm("CEID Growth Factor", "", "1");
        addParm("CEID Max area", "", "100");
        addParm("CEID MF reorientation rate", "", "-1");
        addParm("EpLrcInitial Wall EK", "", "-1");
        addParm("EpLrcInitial Wall CK", "", "-1");
        addParm("EpLrcInitial Shear EK", "", "-1");
        addParm("EpLrcInitial Shear CK", "", "-1");
        addParm("EpLrcInitial Turgor Pressure", "", "-1");
        addParm("EpLrcInitial Growth Factor", "", "1");
        addParm("EpLrcInitial Max area", "", "50");
        addParm("EpLrcInitial MF reorientation rate", "", "-1");
        addParm("LRC Wall EK", "", "-1");
        addParm("LRC Wall CK", "", "-1");
        addParm("LRC Shear EK", "", "-1");
        addParm("LRC Shear CK", "", "-1");
        addParm("LRC Turgor Pressure", "", "-1");
        addParm("LRC Growth Factor", "", "1");
        addParm("LRC Max area", "", "75");
        addParm("LRC MF reorientation rate", "", "-1");
        addParm("Substrate Wall EK", "", "-1");
        addParm("Substrate Wall CK", "", "-1");
        addParm("Substrate Shear EK", "", "-1");
        addParm("Substrate Shear CK", "", "-1");
        addParm("Substrate Turgor Pressure", "", "-1");
        addParm("Substrate Growth Factor", "", "0");
        addParm("Substrate MF reorientation rate", "", "0");
        addParm("Substrate Max area", "", "1000");
        addParm("Source Wall EK", "", "1");
        addParm("Source Wall CK", "", "1");
        addParm("Source Shear EK", "", "1");
        addParm("Source Shear CK", "", "1");
        addParm("Source Turgor Pressure", "", "-1");
        addParm("Source Growth Factor", "", "0");
        addParm("Source MF reorientation rate", "", "0");
        addParm("Source Max area", "", "1000");

    }
    bool initialize(QWidget* parent);
    bool step();
    Mechanics* mechanicsProcess = 0;
    MechanicalGrowth * mechanicalGrowthProcess = 0;
    RootDivide*  divideProcess = 0;
};

class PrintVertexAttr : public Process {
public:
    PrintVertexAttr(const Process& process)
        : Process(process) {
        setName("Model/Root/24 PrintVertexAttr");
        setDesc("PrintVertexAttr.");

    }
    bool step() {
        Mesh* mesh = getMesh("Mesh 1");
        CCIndexDataAttr& indexAttr = mesh->indexAttr();

        Tissue::VertexDataAttr& vMAttr =
            mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");
        CCStructure& cs = mesh->ccStructure(mesh->ccName());
        for(CCIndex v : cs.vertices())
            if(indexAttr[v].selected) {
                mdxInfo << "Vertex: " << v << ", " << indexAttr[v].pos << " - ";
                if(((Tissue::CellData*)vMAttr[v].dualCell != 0))
                    mdxInfo << ((Tissue::CellData*)vMAttr[v].dualCell)->label << endl;
                else
                    mdxInfo << "No dual cell" << endl;
            }
        return false;
    }
};


class HighlightCell : public Process {
public:
    HighlightCell(const Process& process)
        : Process(process) {
        setName("Model/Root/25 Highlight Cell");
        setDesc("HighlightCell.");
        addParm("Face Label", "", "0");

    }
    bool step();
};


class AddFace : public Process {
public:
    AddFace(const Process& process)
        : Process(process) {
        setName("Model/Root/Add Face");
        addParm("Orientation",
                "Orientation",
                "Counter clock-wise",
                QStringList() << "Counter clock-wise"
                              << "Clock-wise");
        addParm("Label", "Label", "0");
    }



    void addFace(CCStructure &cs, CCIndexDataAttr& indexAttr, std::set<CCIndex> vs, int label, int orientation = 0) {

        // the selected orientation
        if(orientation == 0){
            if(parm("Orientation") == "Counter clock-wise")
                orientation = -1;
            else if(parm("Orientation") == "Clock-wise")
                orientation = 1;
            else
                throw(QString("Root::AddFace bad orientation"));
        }
        if(orientation != -1 && orientation != 1)
            throw(QString("Root::AddFace bad orientation"));


        if(vs.size() < 3)
            throw(QString("Root::AddFace needs at least 3 vertices to make a face"));

        // obtain the edges connecting the selected vertices
        std::vector<std::pair<Point3d, Point3d>> polygonSegs;
        for(CCIndex v1 : vs)
            for(CCIndex v2 : vs)
                if(!edgeBetween(cs, v1, v2).isPseudocell())
                    if(std::find(polygonSegs.begin(), polygonSegs.end(), std::make_pair(indexAttr[v1].pos, indexAttr[v2].pos)) == polygonSegs.end() &&
                       std::find(polygonSegs.begin(), polygonSegs.end(), std::make_pair(indexAttr[v2].pos, indexAttr[v1].pos)) == polygonSegs.end())
                            polygonSegs.push_back(std::make_pair(indexAttr[v1].pos, indexAttr[v2].pos));

        // get the centroid
        Point3d centroid;
        int num_v = 0;
        for(CCIndex v : vs) {
            centroid += indexAttr[v].pos;
            num_v++;
        }
        centroid /= num_v;

        // sort the edges and get the sum of the normals, reverse the vertices according to user selection
        // (the polygon can be concave, so we get a majority vote or normals among the edges)
        std::vector<Point3d> ps_ordered = orderPolygonSegs(polygonSegs);
        int nrml_count = 0;
        for(uint i = 1; i < ps_ordered.size(); i++) {
            Point3d p1 = ps_ordered[i-1];
            Point3d p2 = ps_ordered[i];
            Point3d nrml = (p2 - centroid).cross(p1 - centroid);
            nrml_count += (nrml.z() > 0) ? 1 : -1;
        }
        if((nrml_count > 0 && orientation == -1) ||
           (nrml_count < 0 && orientation == 1))
            std::reverse(std::begin(ps_ordered), std::end(ps_ordered));

        // retrieve the original vertices now sorted
        std::vector<CCIndex> vs_final;
        for(Point3d p : ps_ordered)
            for(CCIndex v : vs)
                if(norm(indexAttr[v].pos - p) < EPS)
                    vs_final.push_back(v);

        // add final face
        CCIndex ft = CCIndexFactory.getIndex();
        mdx::addFace(cs, ft, vs_final);
        indexAttr[ft].label = label;

    }


    bool step() {
        Mesh* mesh = getMesh("Mesh 1");
        if(!mesh or mesh->file().isEmpty())
            throw(QString("Root::AddFace No current mesh"));

        QString ccName = mesh->ccName();
        if(ccName.isEmpty())
            throw(QString("Root::AddFace Error, no cell complex selected"));

        CCStructure& cs = mesh->ccStructure(ccName);
        CCIndexDataAttr& indexAttr = mesh->indexAttr();


        // find the selected vertices
        std::set<CCIndex> vs;
        for(CCIndex v : cs.vertices())
            if(indexAttr[v].selected)
                vs.insert(v);

        // the selected orientation
        int orientation = 0;
        if(parm("Orientation") == "Counter clock-wise")
            orientation = -1;
        else if(parm("Orientation") == "Clock-wise")
            orientation = 1;
        else
            throw(QString("Root::AddFace bad orientation"));

        this->addFace(cs, indexAttr, vs, parm("Label").toInt(), orientation);
        mesh->updateAll();

        return false;
    }
};


class DeleteEdges : public Process {
public:
    DeleteEdges(const Process& process)
        : Process(process) {
        setName("Model/Root/Delete Edges");

    }

    bool step() {
        Mesh* mesh = getMesh("Mesh 1");
        if(!mesh or mesh->file().isEmpty())
            throw(QString("Root::AddFace No current mesh"));

        QString ccName = mesh->ccName();
        if(ccName.isEmpty())
            throw(QString("Root::AddFace Error, no cell complex selected"));

        CCStructure& cs = mesh->ccStructure(ccName);
        CCIndexDataAttr& indexAttr = mesh->indexAttr();


        // find the selected vertices
        std::set<CCIndex> vs;
        for(CCIndex v : cs.vertices())
            if(indexAttr[v].selected)
                vs.insert(v);

        for(CCIndex v : vs)
            for(CCIndex vn : cs.neighbors(v))
                if(vs.find(vn) != vs.end()) {
                    CCSignedIndex oe = cs.orientedEdge(v, vn);
                    if(!(~oe).isPseudocell()) {
                        cs.deleteCell(~oe);
                        mdxInfo << "Deleted: " << oe << endl;
                    }
                }
        mesh->updateAll();

        return false;
    }
};



class CreateEdge : public Process {
public:
    CreateEdge(const Process& process)
        : Process(process) {
        setName("Model/Root/Create Edge");

    }

    bool step() {
        Mesh* mesh = getMesh("Mesh 1");
        if(!mesh or mesh->file().isEmpty())
            throw(QString("Root::AddFace No current mesh"));

        QString ccName = mesh->ccName();
        if(ccName.isEmpty())
            throw(QString("Root::AddFace Error, no cell complex selected"));

        CCStructure& cs = mesh->ccStructure(ccName);
        CCIndexDataAttr& indexAttr = mesh->indexAttr();


        // find the selected vertices
        std::vector<CCIndex> vs;
        for(CCIndex v : cs.vertices())
            if(indexAttr[v].selected)
                vs.push_back(v);

        if(vs.size() != 2)
            throw(QString("Too many edges selected"));

        CCIndex e = CCIndexFactory.getIndex();

        cs.addCell(e, +vs[0] -vs[1]);

        mesh->updateAll();

        return false;
    }
};




class ReverseCell : public Process {
public:
    ReverseCell(const Process& process)
        : Process(process) {
        setName("Model/Root/Reverse Cell");

    }

    bool step() {
        Mesh* mesh = getMesh("Mesh 1");
        if(!mesh or mesh->file().isEmpty())
            throw(QString("Root::AddFace No current mesh"));

        QString ccName = mesh->ccName();
        if(ccName.isEmpty())
            throw(QString("Root::AddFace Error, no cell complex selected"));

        CCStructure& cs = mesh->ccStructure(ccName);
        CCIndexDataAttr& indexAttr = mesh->indexAttr();

        for(CCIndex f : cs.faces())
            if(indexAttr[f].selected)
                cs.reverseOrientation(f);

        mesh->updateAll();

        return false;
    }
};




/*
class ReadRootVerticesFromFile : public Process {
public:
    ReadRootVerticesFromFile(const Process &process) : Process(process) {
        setName("Model/Root/Read Root Vertices From File");
        addParm("File", "File", "");
        addParm("Scale", "Scale", "1");
    }
    bool step() {
        Mesh *mesh = getMesh("Mesh 1");
        if (!mesh or mesh->file().isEmpty())
            throw(QString("Root::ReadRootVerticesFromFile No current mesh"));

        CCStructure &cs = mesh->ccStructure("Root");
        CCIndexDataAttr &indexAttr = mesh->indexAttr();

        std::ifstream file(parm("File").toStdString());
        double scale = parm("Scale").toDouble();
        std::string   line;

        while(std::getline(file, line))
        {
            double x = std::stod(line.substr(0, line.find(',')));
            double y = std::stod(line.substr(line.find(',')+1, line.size()));
            CCIndex v = CCIndexFactory.getIndex();
            cs.addCell(v);
            indexAttr[v].pos = Point3d(x/scale, y/scale, 0);
        }
        return false;
    }
};*/

} // namespace ROOT


#endif // ROOT_HPP
