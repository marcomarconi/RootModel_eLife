#ifndef TISSUE_HPP
#define TISSUE_HPP

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <Process.hpp>
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

using namespace mdx;


const double EPS = 1e-6;
const double BIG_VAL = 100000;

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() ;

void SVDDecompX(const Matrix3d &M, Matrix3d &U, Matrix3d &S, Matrix3d &V);

void PolarDecompX(Matrix3d &M,Matrix3d &S, Matrix3d &U);

Point3d Rotate(Point3d v, double angle);

bool PointInPolygon(Point3d point, std::vector<Point3d> points) ;

Point3d shortenLength(Point3d A, float reductionLength) ;

float DistancePtLine(Point3d a, Point3d b, Point3d p) ;

Point3d findClosestLineToLine(Point3d targetLine,
                           Point3d line1, Point3d line2) ;

std::vector<double> softmax(std::vector<double> v) ;

void MBR(std::vector<Point3d> points, std::vector<Point3d>& rect) ;

// extended version of the official one
void neighborhood2D(Mesh& mesh,
                    const CCStructure& cs,
                    const CCIndexDataAttr& indexAttr,
                    std::map<IntIntPair, double>& wallAreas,
                    std::map<IntIntPair, std::set<CCIndex>>& wallEdges);

class Tissue : public Process, public virtual TissueParms {
public:
    Tissue(const Process& process)
        : Process(process) {
        setName("Model/Root/03 Cell Tissue");
        addParm("Set Global Attr Process",
                "Name of Set Global Attributes Process",
                "Model/Root/23 Set Global Attr");
        setIcon(QIcon(":/images/CellDivide.png"));
        addParm("Draw Velocity Rescale", "Draw Velocity Rescale", "0.0");
        addParm("Draw Strain Tensors Rescale", "Draw Strain Tensors Rescale", "0");
        addParm("Draw MF",
                "Draw MF",
                "True",
                QStringList() << "False"
                              << "True");
        addParm("Draw Bounding Box",
                "Draw Bounding Box",
                "False",
                QStringList() << "False"
                              << "True");
        addParm("Draw Auxin-Flux Rescale",
                "Draw Auxin-Flux Rescale",
                 "1");
        addParm("Draw PGD Rescale",
                "Draw PGD Rescale",
                 "1");
        addParm("Draw PINs",
                "Draw PINs",
                 "0.1");
        addParm("Verbose", "Verbose", "False",
                                                        QStringList() << "False"
                                                                      << "True");
    }

    struct CellData;
    struct FaceData;
    struct EdgeData;
    struct VertexData;
    typedef AttrMap<int, CellData> CellDataAttr;
    typedef AttrMap<CCIndex, FaceData> FaceDataAttr;
    typedef AttrMap<CCIndex, EdgeData> EdgeDataAttr;
    typedef AttrMap<CCIndex, VertexData> VertexDataAttr;

    struct Data {
        int dim = -1; // for the chemical solver
    };

    enum VertexType { NoVertexType, Inside, Border, Visual };

    static const char* ToString(VertexType v) {
        switch(v) {
        case NoVertexType:
            return "NoVertexType";
        case Inside:
            return "Inside";
        case Border:
            return "Border";
        case Visual:
            return "Visual";
        default:
            throw(QString("Bad vertex type %1").arg(v));
        }
    }

    // Vertex Data
    struct VertexData : Data {
        // attributes worth saving
        VertexType type = NoVertexType;
        double invmass = -1;
        Point3d restPos,  prevPos,  lastPos;
        Point3d velocity, prevVelocity;
        Data* dualCell = 0; // in case this is a vertex of the dual graph
        bool divisionPoint = false;
        std::map<std::pair<CCIndex, CCIndex>, double> angle;
        // other dynamically loaded attributes
        Point3u dirichlet = Point3u(0, 0, 0); // fixed Dirichlet boundary in x, y, z
        int clusters = 0;
        std::map<const char*, Point3d> corrections;
        Point3d dampedVelocity;
        bool substrate = false, source = false;
        Point3d force;
        // just for debug, these variable have no effect
        std::vector<std::tuple<int, const char*, Point3d>> forces;
        CCIndex V_v0, V_v1, V_e;

        VertexData() {
            dim = 0;
        }

        void resetChems() {}
        void resetMechanics() {
            restPos = prevPos = lastPos; // ?
        }

        void restore(CCIndex v, const CCStructure& cs,
                     const CCIndexDataAttr &indexAttr,const CellDataAttr &cellAttr,
                     const EdgeDataAttr &edgeAttr) {
            type = Tissue::Inside;
            for(CCIndex e : cs.incidentCells(v, 1))
                if(edgeAttr[e].type == Wall)
                    type = Border;
            substrate = false, source = false;
            std::set<int> labels;
            for(CCIndex f : cs.incidentCells(v, 2)) {
                labels.insert(indexAttr[f].label);
                if(cellAttr[indexAttr[f].label].type==Substrate)
                    substrate = true;
                if(cellAttr[indexAttr[f].label].type==Source)
                    source = true;
            }
            clusters = labels.size();
            if(norm(restPos) == 0)
                restPos = indexAttr[v].pos;
        }

        // careful when you change these
        bool read(const QByteArray& ba, size_t& pos) {
            readPOD<VertexType>(type, ba, pos);
            readPOD<Data*>(dualCell, ba, pos);
            readPOD<double>(invmass, ba, pos);
            readPOD<Point3d>(restPos, ba, pos);
            readPOD<Point3d>(prevPos, ba, pos);
            readPOD<Point3d>(lastPos, ba, pos);
            readPOD<Point3d>(velocity, ba, pos);
            readPOD<Point3d>(prevVelocity, ba, pos);
            readPOD<bool>(divisionPoint, ba, pos);
            return true;
        }

        // attribute must be write and read in the same order
        bool write(QByteArray& ba) {
            writePOD<VertexType>(type, ba);
            writePOD<Data*>(dualCell, ba);
            writePOD<double>(invmass, ba);
            writePOD<Point3d>(restPos, ba);
            writePOD<Point3d>(prevPos, ba);
            writePOD<Point3d>(lastPos, ba);
            writePOD<Point3d>(velocity, ba);
            writePOD<Point3d>(prevVelocity, ba);
            writePOD<bool>(divisionPoint, ba);
            return true;
        }

        bool operator==(const VertexData& other) const {
            if(dualCell == other.dualCell)
                return true;
            return false;
        }
    };

    // Cell Data
    enum CellType {
        Undefined,
        QC,
        Columella,
        ColumellaInitial,
        CEI,
        CEID,
        Cortex,
        Endodermis,
        Pericycle,
        Vascular,
        VascularInitial,
        EpLrcInitial,
        Epidermis,
        LRC,
        Substrate,
        Source
    };

    static CellType stringToCellType(const QString& str) {
        if(str == "Undefined")
            return (Undefined);
        else if(str == "QC")
            return (QC);
        else if(str == "Columella")
            return (Columella);
        else if(str == "ColumellaInitial")
            return (ColumellaInitial);
        else if(str == "CEI")
            return (CEI);
        else if(str == "CEID")
            return (CEID);
        else if(str == "Cortex")
            return (Cortex);
        else if(str == "Endodermis")
            return (Endodermis);
        else if(str == "Vascular")
            return (Vascular);
        else if(str == "VascularInitial")
            return (VascularInitial);
        else if(str == "Pericycle")
            return (Pericycle);
        else if(str == "EpLrcInitial")
            return (EpLrcInitial);
        else if(str == "Epidermis")
            return (Epidermis);
        else if(str == "LRC")
            return (LRC);
        else if(str == "Substrate")
            return (Substrate);
        else if(str == "Source")
            return (Source);
        else
            throw(QString("Bad cell type %1").arg(str));
    }

    static const char* ToString(CellType v) {
        switch(v) {
        case Undefined:
            return "Undefined";
        case QC:
            return "QC";
        case Columella:
            return "Columella";
        case ColumellaInitial:
            return "ColumellaInitial";
        case CEI:
            return "CEI";
        case CEID:
            return "CEID";
        case Cortex:
            return "Cortex";
        case Endodermis:
            return "Endodermis";
        case VascularInitial:
            return "VascularInitial";
        case Vascular:
            return "Vascular";
        case Pericycle:
            return "Pericycle";
        case EpLrcInitial:
            return "EpLrcInitial";
        case Epidermis:
            return "Epidermis";
        case LRC:
            return "LRC";
        case Substrate:
            return "Substrate";
        case Source:
            return "Source";
        default:
            throw(QString("Bad cell type %1").arg(v));
        }
    }

    struct CellData : Data {
        // attributes worth saving
        Tissue *tissue = 0;
        int label = -1;
        CellType type = Undefined;
        int lineage = 0;
        CCIndex dualVertex;
        std::set<CCIndex>* cellFaces = 0;
        std::vector<CCIndex> cellVertices;
        std::vector<CCIndex> perimeterFaces, perimeterVertices, perimeterEdges;
        double area = 0, prevArea = 0, restArea = 0;
        double invmassVertices = 1;
        Point3d a1 = Point3d(0, EPS, 0), a2 = Point3d(EPS, 0, 0);
        double mfRORate = -1;
        double growthFactor = 0;
        double lastDivision = 0;
        double divisionCount = 0;
        Point3d axisMin, axisMax, prev_axisMin, prev_axisMax;
        bool periclinalDivision = false;
        int divAlg = -1;
        bool shapeInit = false;
        double cellMaxArea = 1000;
        Point3d restCm;
        Matrix3d invRestMat;
        std::vector<Point3d> restX0;
        double auxin = 0, prevAuxin = 0;
        double Pin1 = 0;
        double Aux1 = 0;
        double PINOID = 0;
        double PP2A = 0;
        double divInhibitor = 0, divPromoter = 0;
        double lifeTime = 0;
        double pressure = 1, pressureMax = -1;
        double auxinProdRate = -1;
        double auxinDecayRate = 0;
        double pinProdRate = -1;
        double pinInducedRate = -1;
        double aux1ProdRate = -1;
        double aux1InducedRate = -1;
        double aux1MaxEdge = -1;

        // other dinamically loaded attribute
        bool selected = false;
        std::pair<int, int> daughters;
        Point3d divVector = Point3d(1., 0., 0.);
        bool mfDelete = false;
        double perimeter = 0;
        double wallStress = -1;
        std::vector<Point3d> box;
        Point3d centroid;
        std::map<CCIndex, double> perimeterAngles;
        double growthRate = 0, axisMin_growthRate = 0, axisMax_growthRate = 0;
        std::vector<double> growthRates, axisMin_grs, axisMax_grs;
        Matrix3d G, E, U, M, M0, S, F = Matrix3d().identity(), R = Matrix3d().identity();
        double gMax = 0, gMin = 0, gAngle = 0, sMax = 0, sMin = 0, sAngle = 0;
        bool divisionAllowed = true;
        double divProb = 0;
        std::map<int, Point3d> auxinFluxes;
        Point3d auxinFluxVector;
        Point3d auxinGradientVector;

        // visuals
        CCIndex a1_v1, a1_v2, a1_e, a2_v1, a2_v2, a2_e,
                auxinFlux_v1, auxinFlux_v2, auxinFlux_v3, auxinFlux_v4, auxinFlux_e, auxinFlux_el, auxinFlux_er, auxinFlux_eb, auxinFlux_f;
        CCIndex axisMin_v, axisMax_v, axisCenter_v, axisMin_e, axisMax_e;
        CCIndex PDGmax_e, PDGmax_v1, PDGmax_v2, PDGmin_e, PDGmin_v1, PDGmin_v2;
        CCIndex box_v[4], box_e[4];

        CellData()
            : cellFaces(new std::set<CCIndex>()) {
            dim = 5;
            box.reserve(4);
            if(type == CEID)
                periclinalDivision = true;
        }

        ~CellData() {

        }
        // BE CAREFUL, every time you change these two functions, you will have to manually reset the attibutes
        bool read(const QByteArray& ba, size_t& pos) {
            readPOD<Tissue*>(tissue, ba, pos);
            readPOD<int>(label, ba, pos);
            readPOD<CellType>(type, ba, pos);
            readPOD<CCIndex>(dualVertex, ba, pos);
            readPOD<double>(area, ba, pos);
            readPOD<double>(prevArea, ba, pos);
            readPOD<double>(restArea, ba, pos);
            readPOD<double>(invmassVertices, ba, pos);
            readPOD<Point3d>(a1, ba, pos);
            readPOD<Point3d>(a2, ba, pos);
            readPOD<double>(mfRORate, ba, pos);
            readPOD<double>(auxin, ba, pos);
            readPOD<double>(Aux1, ba, pos);
            readPOD<double>(Pin1, ba, pos);
            readPOD<double>(auxinProdRate, ba, pos);
            readPOD<double>(auxinDecayRate, ba, pos);
            // cellFaces
            std::vector<CCIndex> v;
            size_t sz;
            readPOD<size_t>(sz, ba, pos);
            v.resize(sz);
            readChar(v.data(), sz * sizeof(CCIndex), ba, pos);
            cellFaces->clear();
            for(CCIndex i : v)
                cellFaces->insert(i);
            readPOD<double>(lastDivision, ba, pos);
            readPOD<bool>(periclinalDivision, ba, pos);
            readPOD<Point3d>(divVector, ba, pos);
            readPOD<double>(cellMaxArea, ba, pos);
            readPOD<Point3d>(axisMax, ba, pos);
            readPOD<Point3d>(axisMin, ba, pos);
            // cellVertices
            std::vector<CCIndex> cv;
            size_t cv_sz;
            readPOD<size_t>(cv_sz, ba, pos);
            cv.resize(cv_sz);
            readChar(cv.data(), cv_sz * sizeof(CCIndex), ba, pos);
            cellVertices.clear();
            for(CCIndex i : cv)
                cellVertices.push_back(i);
            // perimeterFaces
            std::vector<CCIndex> pf;
            size_t pf_sz;
            readPOD<size_t>(pf_sz, ba, pos);
            pf.resize(pf_sz);
            readChar(pf.data(), pf_sz * sizeof(CCIndex), ba, pos);
            perimeterFaces.clear();
            for(CCIndex i : pf)
                perimeterFaces.push_back(i);
            // perimeterEdges
            std::vector<CCIndex> pe;
            size_t pe_sz;
            readPOD<size_t>(pe_sz, ba, pos);
            pe.resize(pe_sz);
            readChar(pe.data(), pe_sz * sizeof(CCIndex), ba, pos);
            perimeterEdges.clear();
            for(CCIndex i : pe)
                perimeterEdges.push_back(i);
            // perimeterVertices
            std::vector<CCIndex> pv;
            size_t pv_sz;
            readPOD<size_t>(pv_sz, ba, pos);
            pv.resize(pv_sz);
            readChar(pv.data(), pv_sz * sizeof(CCIndex), ba, pos);
            perimeterVertices.clear();
            for(CCIndex i : pv)
                perimeterVertices.push_back(i);
            // invRestMatrix
            readPOD<Point3d>(invRestMat[0], ba, pos);
            readPOD<Point3d>(invRestMat[1], ba, pos);
            readPOD<Point3d>(invRestMat[2], ba, pos);
            // restX0
            std::vector<Point3d> rX0;
            size_t rX0_sz;
            readPOD<size_t>(rX0_sz, ba, pos);
            rX0.resize(rX0_sz);
            readChar(rX0.data(), rX0_sz * sizeof(Point3d), ba, pos);
            restX0.clear();
            for(Point3d i : rX0)
                restX0.push_back(i);
            readPOD<Point3d>(restCm, ba, pos);
            //
            readPOD<double>(lifeTime, ba, pos);
            readPOD<double>(pressure, ba, pos);
            readPOD<double>(pressureMax, ba, pos);
            readPOD<bool>(shapeInit, ba, pos);
            readPOD<int>(divAlg, ba, pos);

            return true;
        }

        // attribute must be write and read in the same order
        bool write(QByteArray& ba) {
            writePOD<Tissue*>(tissue, ba);
            writePOD<int>(label, ba);
            writePOD<CellType>(type, ba);
            writePOD<CCIndex>(dualVertex, ba);
            writePOD<double>(area, ba);
            writePOD<double>(prevArea, ba);
            writePOD<double>(restArea, ba);
            writePOD<double>(invmassVertices, ba);
            writePOD<Point3d>(a1, ba);
            writePOD<Point3d>(a2, ba);
            writePOD<double>(mfRORate, ba);
            writePOD<double>(auxin, ba);
            writePOD<double>(Aux1, ba);
            writePOD<double>(Pin1, ba);
            writePOD<double>(auxinProdRate, ba);
            writePOD<double>(auxinDecayRate, ba);
            // cellFaces
            std::vector<CCIndex> v(cellFaces->begin(), cellFaces->end());
            size_t sz = v.size();
            writePOD<size_t>(sz, ba);
            writeChar(v.data(), sz * sizeof(CCIndex), ba);
            writePOD<double>(lastDivision, ba);
            writePOD<bool>(periclinalDivision, ba);
            writePOD<Point3d>(divVector, ba);
            writePOD<double>(cellMaxArea, ba);
            writePOD<Point3d>(axisMax, ba);
            writePOD<Point3d>(axisMin, ba);
            // cellVertices
            std::vector<CCIndex> cv(cellVertices.begin(), cellVertices.end());
            size_t cv_sz = cv.size();
            writePOD<size_t>(cv_sz, ba);
            writeChar(cv.data(), cv_sz * sizeof(CCIndex), ba);
            // perimeterFaces
            std::vector<CCIndex> pf(perimeterFaces.begin(), perimeterFaces.end());
            size_t pf_sz = pf.size();
            writePOD<size_t>(pf_sz, ba);
            writeChar(pf.data(), pf_sz * sizeof(CCIndex), ba);
            // perimeterEdges
            std::vector<CCIndex> pe(perimeterEdges.begin(), perimeterEdges.end());
            size_t pe_sz = pe.size();
            writePOD<size_t>(pe_sz, ba);
            writeChar(pe.data(), pe_sz * sizeof(CCIndex), ba);
            // perimeterVertices
            std::vector<CCIndex> pv(perimeterVertices.begin(), perimeterVertices.end());
            size_t pv_sz = pv.size();
            writePOD<size_t>(pv_sz, ba);
            writeChar(pv.data(), pv_sz * sizeof(CCIndex), ba);
            // invRestMatrix
            writePOD<Point3d>(invRestMat[0], ba);
            writePOD<Point3d>(invRestMat[1], ba);
            writePOD<Point3d>(invRestMat[2], ba);
            // restX0
            std::vector<Point3d> rx0(restX0.begin(), restX0.end());
            size_t rx0_sz = rx0.size();
            writePOD<size_t>(rx0_sz, ba);
            writeChar(rx0.data(), rx0_sz * sizeof(Point3d), ba);
            writePOD<Point3d>(restCm, ba);
            //
            writePOD<double>(lifeTime, ba);
            writePOD<double>(pressure, ba);
            writePOD<double>(pressureMax, ba);
            writePOD<bool>(shapeInit, ba);
            writePOD<int>(divAlg, ba);

            return true;
        }

        bool operator==(const CellData& other) const {
            if(area == other.area)
                return true;
            return false;
        }


        // reset chemicals
        void resetChems() {
            auxin = 0;
            Pin1 = 0;
            Aux1 = 0;
            PP2A = 0;
            PINOID = 0;
            divPromoter = 0;
            divInhibitor = 0;
            auxinFluxes.clear();
            auxinFluxVector = Point3d(0, 0, 0);
        }

        void resetMechanics(FaceDataAttr& faceAttr) {
            a1 = Point3d(1, 1, 0) * EPS, a2 = Point3d(1, 1, 0) * EPS;
            restArea = area;
            restCm = centroid;
            pressureMax = -1;
            shapeInit = false;
            invmassVertices = 1;            
            mfRORate = -1;
            wallStress = -1;
            for(CCIndex f : *cellFaces){
                FaceData &fD = faceAttr[f];
                fD.resetMechanics();
            }
            growthRates.clear();
        }

        void updateGeom(const CCIndexDataAttr& indexAttr, FaceDataAttr& faceAttr, EdgeDataAttr& edgeAttr, double Dt) {

            prevArea = area;
            area = 0;
            centroid = Point3d(0, 0, 0);
            // deformation gradient and green strain tensor
            G = 0, E = 0, F = 0;
            int i = 0;
            for(CCIndex f : *cellFaces) {
                FaceData fD = faceAttr[f];
                area += indexAttr[f].measure;
                centroid += indexAttr[f].pos;
                if(fD.type == Membrane) {
                    E += fD.E;
                    G += fD.G;
                    F += fD.F;
                    i++;
                }
            }
            centroid /= cellFaces->size();
            E /= i;
            G /= i;
            F /= i;

            // perimeter length
            perimeter = 0;
            for(CCIndex e : perimeterEdges)
                perimeter += edgeAttr[e].length;

            // cell growth rate
            if(Dt > 0)
                growthRates.push_back((area - prevArea) / (prevArea * Dt));
            if(growthRates.size() > 100)
                growthRates.erase(growthRates.begin());
            growthRate = 0;
            int grs = 1;
            for(double gr : growthRates)
                if(gr > 0) {
                    growthRate += gr;
                    grs++;
                }
            growthRate /= grs;

            // principal components of the green strain tensor
            gAngle = gMax = gMin = 0;
            if((G[0][0] - G[1][1]) != 0) {
                gAngle = (0.5 * atan(2 * G[0][1] / (G[0][0] - G[1][1]) )) ;
                double left = (G[0][0] + G[1][1]) / 2;
                double right = sqrt(pow((G[0][0] - G[1][1]) / 2, 2) + pow(G[0][1]/2, 2));
                gMax = left + right;
                gMin = left - right;
            }
            if(G[0][0] < G[1][1])
                if(gAngle > -(M_PI/4) && gAngle < (M_PI/4))
                    gAngle += (M_PI/2);


        }

        void restore(Mesh* mesh, Tissue* tissue,
                     const CCStructure& cs,
                     const CCStructure& csDual,
                     CCIndexDataAttr& indexAttr,
                     FaceDataAttr& faceAttr,
                     EdgeDataAttr& edgeAttr,
                     VertexDataAttr& vMAttr) {
            this->tissue = tissue;
            if(lineage == 0)
                lineage = label;
            // find my faces
            uint old_cellFaces = cellFaces->size();
            cellFaces->clear();
            for(CCIndex f : cs.faces())
                if(indexAttr[f].label == label) {
                    cellFaces->insert(f);
                    faceAttr[f].owner = this;
                }
            if(old_cellFaces != cellFaces->size())
                     shapeInit = false;

            // update the geometry
            updateGeom(indexAttr, faceAttr, edgeAttr, 0);

            // find my dual vertex
            if(csDual.vertices().size() > 0) {
                //throw(QString("CellData::restore empty dual graph"));
                CCIndex closest = csDual.vertices()[0];
                for(CCIndex v : csDual.vertices())
                        if(norm(indexAttr[v].pos - centroid) < norm(indexAttr[closest].pos - centroid))
                            closest = v;
                indexAttr[closest].label = label;
                vMAttr[closest].dualCell = this;
                dualVertex = closest;
            }


            // set perimeter edges and faces, and set average restLength for edges
            cellVertices.clear();
            perimeterFaces.clear();
            perimeterEdges.clear();
            std::vector<CCIndex> tmpEdges, tmpVertices;
            perimeterVertices.clear();
            for(CCIndex f : *cellFaces)
                for(CCIndex e : cs.incidentCells(f, 1)) {
                    if(std::find(cellVertices.begin(), cellVertices.end(), cs.edgeBounds(e).first) == cellVertices.end())
                        cellVertices.push_back(cs.edgeBounds(e).first);
                    if(std::find(cellVertices.begin(), cellVertices.end(), cs.edgeBounds(e).second) == cellVertices.end())
                        cellVertices.push_back(cs.edgeBounds(e).second);
                    if(vMAttr[cs.edgeBounds(e).first].invmass == -1)
                        vMAttr[cs.edgeBounds(e).first].invmass = invmassVertices;
                    if(vMAttr[cs.edgeBounds(e).second].invmass == -1)
                        vMAttr[cs.edgeBounds(e).second].invmass = invmassVertices;
                    std::set<CCIndex> edgeFaces = cs.incidentCells(e, 2);
                    CCIndex f1 = *(edgeFaces.begin());
                    CCIndex f2;
                    if(edgeFaces.size() > 1)
                        f2 = *(++edgeFaces.begin());
                    if(mesh->getLabel(indexAttr[f1].label) == mesh->getLabel(indexAttr[f2].label))
                        continue;
                    if(std::find(perimeterFaces.begin(), perimeterFaces.end(), f) == perimeterFaces.end())
                        perimeterFaces.push_back(f);
                    if(std::find(tmpEdges.begin(), tmpEdges.end(), e) == tmpEdges.end())
                        tmpEdges.push_back(e);
                    if(std::find(tmpVertices.begin(), tmpVertices.end(), cs.edgeBounds(e).first) == tmpVertices.end())
                        tmpVertices.push_back(cs.edgeBounds(e).first);
                    if(std::find(tmpVertices.begin(), tmpVertices.end(), cs.edgeBounds(e).second) == tmpVertices.end())
                        tmpVertices.push_back(cs.edgeBounds(e).second);
                }

            // sort the perimeter edges and vertices
            std::vector<std::pair<Point3d,Point3d> > polygonSegs;
            std::vector<Point3d>  polygonPoints;
            for(CCIndex e : tmpEdges) {
                auto eb = cs.edgeBounds(e);
                polygonSegs.push_back(make_pair(indexAttr[eb.first].pos, indexAttr[eb.second].pos));
            }
            for(CCIndex v : tmpVertices)
                polygonPoints.push_back(indexAttr[v].pos);
            std::vector<std::pair<Point3d,Point3d> > polygonSegsOrig = polygonSegs;
            std::vector<Point3d>  polygonPointsOrdered = mdx::orderPolygonSegs(polygonSegs, true);
            for(auto p : polygonSegs) {
                auto it = std::find(polygonSegsOrig.begin(), polygonSegsOrig.end(), p);
                if(it == polygonSegsOrig.end()) {
                    Point3d tmp = p.first;
                    p.first = p.second;
                    p.second = tmp;
                    it = std::find(polygonSegsOrig.begin(), polygonSegsOrig.end(), p);
                }
                perimeterEdges.push_back(tmpEdges[std::distance(polygonSegsOrig.begin(), std::find(polygonSegsOrig.begin(), polygonSegsOrig.end(), p))]);
            }
            for(auto p : polygonPointsOrdered)
                perimeterVertices.push_back(tmpVertices[std::distance(polygonPoints.begin(), std::find(polygonPoints.begin(), polygonPoints.end(), p))]);

            // set the perimeter vertices' angles
            int n = perimeterVertices.size();
            for(int i = 0; i < n; i++) {
                CCIndex v = perimeterVertices[i];
                if(perimeterAngles.find(v) == perimeterAngles.end()){
                    int prev = i-1;
                    if(prev == -1)
                        prev = n-1;
                    int next = i+1;
                    if(next == n)
                        next = 0;
                    Point3d p0 = indexAttr[perimeterVertices[prev]].pos - indexAttr[perimeterVertices[i]].pos;
                    Point3d p1 = indexAttr[perimeterVertices[next]].pos - indexAttr[perimeterVertices[i]].pos;
                    perimeterAngles[v] = mdx::angle(p0, p1);
                }
                CCIndex e1, e2;
                for(CCIndex e : cs.incidentCells(v, 1))
                    if(find(perimeterEdges.begin(), perimeterEdges.end(), e) != perimeterEdges.end()){
                        if(e1.isPseudocell())
                            e1 = e;
                        else
                            e2 = e;
                    }
                if(e1.isPseudocell() || e2.isPseudocell())
                    throw(QString("Something wrong here"));
                int prev = i-1;
                if(prev == -1) prev = n-1;
                int next = i+1;
                if(next == n)  next = 0;
                Point3d p0 = indexAttr[perimeterVertices[prev]].pos - indexAttr[perimeterVertices[i]].pos;
                Point3d p1 = indexAttr[perimeterVertices[next]].pos - indexAttr[perimeterVertices[i]].pos;
                if(vMAttr[v].angle[make_pair(e1, e2)] == 0)
                    vMAttr[v].angle[make_pair(e1, e2)] = mdx::angle(p0, p1);
            }
        }

        void
        update(const CCStructure& cs, const CCIndexDataAttr& indexAttr, FaceDataAttr& faceAttr, EdgeDataAttr &edgeAttr, VertexDataAttr& vMAttr, double Dt) {

            // Geometry
            updateGeom(indexAttr, faceAttr, edgeAttr, Dt);

            // Misc updates
            lastDivision += Dt;

            if(perimeterVertices.size() > 0) {
                // bounding box
                std::vector<Point3d> points;
                for(auto v : perimeterVertices)
                    points.push_back(indexAttr[v].pos);
                MBR(points, box);
                prev_axisMax = axisMax;
                prev_axisMin = axisMin;
                axisMax = box[0] - box[1];
                axisMin = box[1] - box[2];
                if(norm(axisMax) < norm(axisMin)) {
                    Point3d tmp = axisMin;
                    axisMin = axisMax;
                    axisMax = tmp;
                }
                // Axis growth rates, needed for PGD visualization
                if(Dt > EPS) {
                    double axisMax_gr = (norm(axisMax) - norm(prev_axisMax)) / Dt;
                    double axisMin_gr = (norm(axisMin) - norm(prev_axisMin)) / Dt;
                    if(axisMax_gr > 0 && axisMin_gr > 0) {
                        axisMax_grs.insert(axisMax_grs.begin(), axisMax_gr);
                        axisMin_grs.insert(axisMin_grs.begin(), axisMin_gr);
                    }
                    uint max_values = 30;
                    if(axisMax_grs.size() > max_values || axisMin_grs.size() > max_values) {
                        axisMax_grs.resize(max_values);
                        axisMin_grs.resize(max_values);
                    }
                    axisMax_growthRate = 0, axisMin_growthRate = 0;
                    for(uint i = 0; i < axisMax_grs.size(); i++) {
                        axisMax_growthRate += axisMax_grs[i];
                        axisMin_growthRate += axisMin_grs[i] ;
                    }
                    if(axisMax_grs.size() > 0)
                        axisMax_growthRate /= axisMax_grs.size();
                    if(axisMin_grs.size() > 0)
                        axisMin_growthRate /= axisMin_grs.size();
                }
                /*
                // rest shape tensor
                M0 = 0;
                for(auto v : perimeterVertices) {
                    M0[0][0] += pow(vMAttr[v].restPos[0] - restCm[0], 2);
                    M0[1][1] += pow(vMAttr[v].restPos[1] - restCm[1], 2);
                    M0[0][1] += (vMAttr[v].restPos[0] - restCm[0]) * (vMAttr[v].restPos[1] - restCm[1]);
                    M0[1][0] += (vMAttr[v].restPos[0] - restCm[0]) * (vMAttr[v].restPos[1] - restCm[1]);
                }
                M0 /= perimeterVertices.size();

                // current shape tensor
                M = 0;
                for(auto v : perimeterVertices) {
                    M[0][0] += pow(indexAttr[v].pos[0] - centroid[0], 2);
                    M[1][1] += pow(indexAttr[v].pos[1] - centroid[1], 2);
                    M[0][1] += (indexAttr[v].pos[0] - centroid[0]) * (indexAttr[v].pos[1] - centroid[1]);
                    M[1][0] += (indexAttr[v].pos[0] - centroid[0]) * (indexAttr[v].pos[1] - centroid[1]);
                }
                M /= perimeterVertices.size();

                // shape tensor strain
                S = (M - M0);
                for(int i = 0; i < 3; i++)
                    for(int j = 0; j < 3; j++)
                        S[i][j] = S[i][j] / transpose(M0)[i][j];

                // principal components of the shape strain tensor
                sAngle = sMax = sMin = 0;
                if((S[0][0] - S[1][1]) != 0) {
                    sAngle = (0.5 * atan(2 * S[0][1] / (S[0][0] - S[1][1]) )) ;
                    double left = (S[0][0] + S[1][1]) / 2;
                    double right = sqrt(pow((S[0][0] - S[1][1]) / 2, 2) + pow(S[0][1]/2, 2));
                    sMax = left + right;
                    sMin = left - right;
                }
                if(S[0][0] < S[1][1])
                    if(sAngle > -(M_PI/4) && sAngle < (M_PI/4))
                        sAngle += (M_PI/2);
                */

                /*if(axisMax[1] < 0)
                    axisMax *= -1;
                if(axisMin[1] < 0)
                    axisMin *= -1;*/
            } else
                mdxInfo << "Cell " << label << " has no vertices!?!?" << endl;

        }

        void setType(CellType newType, VertexDataAttr& vMAttr) {
            type = newType;
            for(CCIndex v : cellVertices) {
                vMAttr[v].substrate = newType == Substrate;
                vMAttr[v].source = newType == Source;
            }

        }

        void division(const CCStructure &cs, CellDataAttr& cellAttr,
                       FaceDataAttr& faceAttr, EdgeDataAttr &edgeAttr, CellData& cD1, CellData& cD2, std::map<CellType, int> maxAreas, bool ignoreCellType=false) ;


    };

    // Edge Data
    enum EdgeType { NoEdgeType, Shear, Wall };

    static const char* ToString(EdgeType v) {
        switch(v) {
        case Shear:
            return "Shear";
        case Wall:
            return "Wall";
        case NoEdgeType:
            return "NoEdgeType";
        default:
            throw(QString("Bad edge type %1").arg(v));
        }
    }

    struct EdgeData : Data {
        // attributes worth saving
        EdgeType type = NoEdgeType;
        double restLength = 0;
        double prevLength = 0, length = 0;
        double prev_strain = 0, strain = 0;
        double strainRate = 0;
        double eStiffness = -1;
        double cStiffness = -1;
        double intercellularAuxin = 0;
        std::map<int, double> Aux1;
        std::map<int, double> Pin1;
        std::map<int, double> PINOID;
        std::map<int, double> PP2A;
        // other dynamical attrs
        std::map<int, double> auxinRatio;
        std::map<int, double> auxinGrad;
        Point3d  midPoint;
        std::map<int, double> angle;
        std::map<CCIndex, Point3d> outwardNormal;
        std::map<int, double> auxinFluxImpact, MFImpact, geomImpact;
        std::map<int, double> pin1Sensitivity, pin1SensitivityRaw;
        double sigmaEv = 0, sigmaEe = 0;
        Point3d sigmaForce = Point3d(0, 0, 0);
        std::vector<Point3d> sigmaForces;
        Point3d pressureForce = Point3d(0, 0, 0);

        CCIndex pin_v1_1, pin_v2_1, pin_v3_1, pin_v4_1, pin_e1_1, pin_e2_1, pin_e3_1, pin_e4_1, pin_f_1;
        CCIndex pin_v1_2, pin_v2_2, pin_v3_2, pin_v4_2, pin_e1_2, pin_e2_2, pin_e3_2, pin_e4_2, pin_f_2;

        EdgeData() {
            dim = 1;
        } // for the chemical solver

        // reset chemicals
        void resetChems() {
            intercellularAuxin = 0;
            Aux1.clear();
            Pin1.clear();
            PINOID.clear();
            PP2A.clear();
            auxinFluxImpact.clear();
            pin1Sensitivity.clear();
        }
        void resetMechanics() {
            restLength = length;
            prevLength = length;
        }

        void restore(CCIndex e,
                     const CCStructure& cs,
                     const CCIndexDataAttr& indexAttr,
                     const CellDataAttr& cellAttr) {
            auto eb = cs.edgeBounds(e);
            length = norm(indexAttr[eb.first].pos - indexAttr[eb.second].pos);
            if(restLength <= 0)
                restLength = length;

            // determine whether we are a wall or internal edge
            std::vector<int> incident_labels;
            for(CCIndex f : cs.incidentCells(e, 2)) {
                int label = indexAttr[f].label;
                incident_labels.push_back(label);
            }
            if((incident_labels.size() == 2 && incident_labels[0] != incident_labels[1]) ||
               cs.onBorder(e))
                this->type = Tissue::Wall;
            else
                this->type = Tissue::Shear;


        }

        // to be called at each step
        void update(CCIndex e,
                    const CCStructure& cs,
                    const CCIndexDataAttr& indexAttr,
                    VertexDataAttr& vMAttr,
                    FaceDataAttr& faceAttr) {
            auto eb = cs.edgeBounds(e);
            Point3d vPos = indexAttr[eb.first].pos;
            Point3d nPos = indexAttr[eb.second].pos;
            prevLength = length;
            length = norm(vPos - nPos);
            prev_strain = strain;
            strain = (length - restLength) / restLength;
            Point3d versor = vPos - nPos;
            Point3d dv = vMAttr[eb.first].velocity - vMAttr[eb.second].velocity;
            strainRate = (versor * dv) / prevLength;

            midPoint = (nPos + (vPos - nPos) / 2.);
            outwardNormal.clear();
            std::vector<int> incident_labels;
            for(CCIndex f : cs.incidentCells(e, 2)) {
                incident_labels.push_back(indexAttr[f].label);
                Point3d versor = (midPoint - indexAttr[f].pos);
                outwardNormal[f] = Point3d(0., 0., 1.) ^ (vPos - nPos);
                outwardNormal[f] *= 1. / norm(outwardNormal[f]);
                if(outwardNormal[f] * versor < 0.)
                    outwardNormal[f] *= -1.;
                Point3d eVector = vPos - nPos;
                // angle is relative to MF a1 orientation
                /*
                 *               90º
                 *                ^
                 *                |
                 *                |
                 *                |
                 *                |
                 * 0º <-----------o---------------> 180º
                 *                |
                 *                |
                 *                |
                 *                |
                 *                v
                 *               90º
                 */

                FaceData& fD = faceAttr[f];
                if(!fD.owner)
                    throw(QString("Tissue::EdgeData::update: orphan face, something wrong with tissue restore or cell division?"));
                // angle[fD.owner->label]  = fabs(atan2(eVector[1], eVector[0])) *
                // 180./M_PI;
                angle[indexAttr[f].label] =
                    acos((eVector * fD.a1) / (norm(eVector) * norm(fD.a1))) *
                    (180. / M_PI);
            }
        }

        bool read(const QByteArray& ba, size_t& pos) {
            readPOD<EdgeType>(type, ba, pos);
            readPOD<double>(length, ba, pos);
            readPOD<double>(prevLength, ba, pos);
            readPOD<double>(prev_strain, ba, pos);
            readPOD<double>(strain, ba, pos);
            readPOD<double>(strainRate, ba, pos);
            readPOD<double>(restLength, ba, pos);
            readPOD<double>(intercellularAuxin, ba, pos);
            readPOD<double>(eStiffness, ba, pos);
            readPOD<double>(cStiffness, ba, pos);
            // PIN
            std::vector<int> k; std::vector<double>v;
            size_t szk, szv;
            readPOD<size_t>(szk, ba, pos);
            readPOD<size_t>(szv, ba, pos);
            k.resize(szk);
            v.resize(szv);
            readChar(k.data(), szk * sizeof(int), ba, pos);
            readChar(v.data(), szv * sizeof(double), ba, pos);
            Pin1.clear();
            for(uint i = 0; i < k.size(); i++)
                Pin1[k[i]] = v[i];
            // AUX1
            k.clear(); v.clear();
            readPOD<size_t>(szk, ba, pos);
            readPOD<size_t>(szv, ba, pos);
            k.resize(szk);
            v.resize(szv);
            readChar(k.data(), szk * sizeof(int), ba, pos);
            readChar(v.data(), szv * sizeof(double), ba, pos);
            Aux1.clear();
            for(uint i = 0; i < k.size(); i++)
                Aux1[k[i]] = v[i];
            // PP2A
            k.clear(); v.clear();
            readPOD<size_t>(szk, ba, pos);
            readPOD<size_t>(szv, ba, pos);
            k.resize(szk);
            v.resize(szv);
            readChar(k.data(), szk * sizeof(int), ba, pos);
            readChar(v.data(), szv * sizeof(double), ba, pos);
            PP2A.clear();
            for(uint i = 0; i < k.size(); i++)
                PP2A[k[i]] = v[i];
            // PINOID
            k.clear(); v.clear();
            readPOD<size_t>(szk, ba, pos);
            readPOD<size_t>(szv, ba, pos);
            k.resize(szk);
            v.resize(szv);
            readChar(k.data(), szk * sizeof(int), ba, pos);
            readChar(v.data(), szv * sizeof(double), ba, pos);
            PINOID.clear();
            for(uint i = 0; i < k.size(); i++)
                PINOID[k[i]] = v[i];
            return true;
        }

        bool write(QByteArray& ba) {
            writePOD<EdgeType>(type, ba);
            writePOD<double>(length, ba);
            writePOD<double>(prevLength, ba);
            writePOD<double>(prev_strain, ba);
            writePOD<double>(strain, ba);
            writePOD<double>(strainRate, ba);
            writePOD<double>(restLength, ba);
            writePOD<double>(intercellularAuxin, ba);
            writePOD<double>(eStiffness, ba);
            writePOD<double>(cStiffness, ba);
            // PIN
            std::vector<int> k; std::vector<double>v;
            for( auto it = Pin1.begin(); it != Pin1.end(); ++it ) {
                    v.push_back( it->second );
                    k.push_back(it->first);
                }
            size_t szv = v.size();
            size_t szk = k.size();
            writePOD<size_t>(szk, ba);
            writePOD<size_t>(szv, ba);
            writeChar(k.data(), szk * sizeof(int), ba);
            writeChar(v.data(), szv * sizeof(double), ba);
            // AUX1
            k.clear(); v.clear();
            for( auto it = Aux1.begin(); it != Aux1.end(); ++it ) {
                    v.push_back( it->second );
                    k.push_back(it->first);
                }
            szv = v.size();
            szk = k.size();
            writePOD<size_t>(szk, ba);
            writePOD<size_t>(szv, ba);
            writeChar(k.data(), szk * sizeof(int), ba);
            writeChar(v.data(), szv * sizeof(double), ba);
            // PINOID
            k.clear(); v.clear();
            for( auto it = PINOID.begin(); it != PINOID.end(); ++it ) {
                    v.push_back( it->second );
                    k.push_back(it->first);
                }
            szv = v.size();
            szk = k.size();
            writePOD<size_t>(szk, ba);
            writePOD<size_t>(szv, ba);
            writeChar(k.data(), szk * sizeof(int), ba);
            writeChar(v.data(), szv * sizeof(double), ba);
            // PP2A
            k.clear(); v.clear();
            for( auto it = PP2A.begin(); it != PP2A.end(); ++it ) {
                    v.push_back( it->second );
                    k.push_back(it->first);
                }
            szv = v.size();
            szk = k.size();
            writePOD<size_t>(szk, ba);
            writePOD<size_t>(szv, ba);
            writeChar(k.data(), szk * sizeof(int), ba);
            writeChar(v.data(), szv * sizeof(double), ba);
            return true;
        }

        bool operator==(const EdgeData& other) const {
            if(restLength == other.restLength and
               cStiffness == other.cStiffness)
                return true;
            return false;
        }
    };

    // FaceData
    enum FaceType { NoFaceType, Internal, Membrane };


    static const char* ToString (FaceType v) {
        switch(v) {
        case Internal:
            return "Internal";
        case Membrane:
            return "Membrane";
        case NoFaceType:
            return "NoFaceType";
        default:
            throw(QString("Bad edge type %1").arg(v));
        }
    }

    struct FaceData : Data {
        // attributes worth saving
        FaceType type = NoFaceType;
        Point3d a1 = Point3d(0, 1, 0), a2 = Point3d(0, -1, 0);
        double restAreaFace = 0;
        Point3d restPos[3];
        Matrix2d invRestMat;
        CellData* owner = 0;
        double area = 0;
        double prevArea = 0;
        // dynamical attributes
        // just a copy of owner's, only relevant for visuals
        double stress = 0;
        double pressure = 0;
        double edgeStrain = 0;
        double edgeStiffness = 0;
        double MFImpact = 0, auxinFluxImpact = 0, geomImpact = 0, auxinRatio = 0, auxinGrad = 0;
        double pin1Sensitivity = 0, pin1SensitivityRaw = 0;
        Matrix3d E, G, V, dV, F = Matrix3d().identity(), R = Matrix3d().identity(), U = Matrix3d().identity(), Ev, F1, F2, sigmaA;
        Point3d divVector = Point3d(1., 0., 0.);
        double auxin = 0, Aux1Cyt = 0, Pin1Cyt = 0, divInhibitorCyt = 0,
               Aux1Mem = 0, Pin1Mem = 0, intercellularAuxin = 0, PINOIDMem = 0, PP2AMem = 0,
               growthRate = 0;
        CCIndex G_v0, G_v1, G_v2, G_e1, G_e2, a1_v1, a1_v2, a1_e, a2_v1, a2_v2, a2_e;

        FaceData() {
            dim = 2;
            E[0] = Point3d(1., 0., 0.);
            E[1] = Point3d(0., 1., 0.);
            E[2] = Point3d(0., 0., 1.);
        }

        // reset chemicals, not really needed as faces just copy owner's
        void resetChems() {}
        void resetMechanics() {
            a1 = Point3d(0, EPS, 0), a2 = Point3d(EPS, 0, 0);
            restAreaFace = area;
            for(Point3d p : restPos)
                p = 0;
            invRestMat = 0;
        }


        void restore(CCIndex f, const CCStructure& cs, const CCIndexDataAttr& indexAttr) {
            a1 = owner->a1;
            a2 = owner->a2;
            area = indexAttr[f].measure;
            if(restAreaFace == 0)
                restAreaFace = area;
            type = Internal;
            if(cs.onBorder(f))
                type = Membrane;
            for(CCIndex fn : cs.neighbors(f))
                if(indexAttr[fn].label != indexAttr[f].label)
                    type = Tissue::Membrane;
        }

        void update(CCIndex f,
                    const CCStructure &cs, const CCIndexDataAttr& indexAttr, VertexDataAttr &vMAttr) {
            prevArea = area;
            area = indexAttr[f].measure;
            growthRate = owner->growthRate;
            std::vector<CCIndex> vs = faceVertices(cs, f);
            Point3d x_p[3] = {vMAttr[vs[0]].prevPos,vMAttr[vs[1]].prevPos,vMAttr[vs[2]].prevPos}; ////// NW
            Point3d X_p[3] = {indexAttr[vs[0]].pos, indexAttr[vs[1]].pos, indexAttr[vs[2]].pos};
            F = DefGradient(x_p, X_p);
            F[2][2] = 1;
            G = 0.5 * (transpose(F) * F - F.identity());
            //Point3d v_p[3] = {vMAttr[vs[0]].prevVelocity,vMAttr[vs[1]].prevVelocity,vMAttr[vs[2]].prevVelocity};
            Point3d V_p[3] = {vMAttr[vs[0]].velocity,vMAttr[vs[1]].velocity,vMAttr[vs[2]].velocity};
            V = DefGradient(X_p, V_p);
            E = 0.5 * (transpose((V) + V));
            //PolarDecompX(F, R, U);
        }

        bool read(const QByteArray& ba, size_t& pos) {
            readPOD<Point3d>(a1, ba, pos);
            readPOD<Point3d>(a2, ba, pos);
            readPOD<double>(restAreaFace, ba, pos);
            readPOD<double>(area, ba, pos);
            readPOD<double>(prevArea, ba, pos);
            readPOD<FaceType>(type, ba, pos);
            readPOD<CellData*>(owner, ba, pos);
            readPOD<Point3d>(restPos[0], ba, pos);
            readPOD<Point3d>(restPos[1], ba, pos);
            readPOD<Point3d>(restPos[2], ba, pos);
            readPOD<Point2d>(invRestMat[0], ba, pos);
            readPOD<Point2d>(invRestMat[1], ba, pos);
            readPOD<Point2d>(invRestMat[2], ba, pos);
            return true;
        }

        bool write(QByteArray& ba) {
            writePOD<Point3d>(a1, ba);
            writePOD<Point3d>(a2, ba);
            writePOD<double>(restAreaFace, ba);
            writePOD<double>(area, ba);
            writePOD<double>(prevArea, ba);
            writePOD<FaceType>(type, ba);
            writePOD<CellData*>(owner, ba);
            writePOD<Point3d>(restPos[0], ba);
            writePOD<Point3d>(restPos[1], ba);
            writePOD<Point3d>(restPos[2], ba);
            writePOD<Point2d>(invRestMat[0], ba);
            writePOD<Point2d>(invRestMat[1], ba);
            writePOD<Point2d>(invRestMat[2], ba);
            return true;
        }
        bool operator==(const FaceData& other) const {
            if(area == other.area)
                return true;
            return false;
        }
    };

    class Subdivide : virtual public mdx::Subdivide {

    public:
        Subdivide() {}
        Subdivide(CCIndexDataAttr& _indexAttr,
                  Tissue::VertexDataAttr& _vtxAttr,
                  Tissue::EdgeDataAttr& _edgeAttr,
                  Tissue::FaceDataAttr& _faceAttr,
                  Tissue::CellDataAttr& _cellAttr)
            : indexAttr(&_indexAttr)
            , cellAttr(&_cellAttr)
            , edgeAttr(&_edgeAttr)
            , faceAttr(&_faceAttr)
            , vtxAttr(&_vtxAttr) {}

        void splitCellUpdate(Dimension dim,
                             const CCStructure& cs,
                             const CCStructure::SplitStruct& ss,
                             CCIndex otherP = CCIndex(),
                             CCIndex otherN = CCIndex(),
                             double interpPos = 0.5, bool cellDivision=false);

    private:
        Mesh* mesh = 0;
        CCIndexDataAttr* indexAttr = 0;
        CellDataAttr* cellAttr = 0;
        EdgeDataAttr* edgeAttr = 0;
        FaceDataAttr* faceAttr = 0;
        VertexDataAttr* vtxAttr = 0;
    };

    // Initialize tissue, called from GUI thread
    bool initialize(QWidget* parent) {
        Mesh* mesh = getMesh("Mesh 1");
        if(!mesh)
            throw(QString("CellTissueProcess::initialize No current mesh"));

        return initialize();
    }
    bool initialize(bool create_dual = true, bool extended_dual = true);
    void restore(CCStructure csCurr);
    bool step(double Dt);
    void createDualExtended(CCStructure &cs, CCStructure &csDual);

    CellTissue& tissue() {
        return cellTissue;
    }

    const QString& tissueName() const {
        return TissueName;
    }
    const QString& tissueDualName() const {
        return TissueDualName;
    }

    void resetChems() {
        Mesh* mesh = getMesh("Mesh 1");
        if(!mesh)
            throw(QString("Tissue::resetChems No current mesh"));
        Tissue::CellDataAttr& cellAttr =
           mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
        Tissue::EdgeDataAttr& edgeAttr =
           mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
        Tissue::FaceDataAttr& faceAttr =
           mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");

        CCStructure& cs =mesh->ccStructure("Tissue");

        for(auto c : cellAttr) {
            Tissue::CellData& cD = cellAttr[c.first];
            cD.resetChems();
        }
        for(CCIndex f : cs.faces()) {
            Tissue::FaceData& fD = faceAttr[f];
            fD.resetChems();
        }
        for(CCIndex e : cs.edges()) {
            Tissue::EdgeData& eD = edgeAttr[e];
            eD.resetChems();
        }
        intercellularChems.clear();
    }

    void resetMechanics() {
        Mesh* mesh = getMesh("Mesh 1");
        if(!mesh)
            throw(QString("Tissue::resetMechanics No current mesh"));
        Tissue::CellDataAttr& cellAttr =
           mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
        Tissue::EdgeDataAttr& edgeAttr =
           mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
        Tissue::FaceDataAttr& faceAttr =
           mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");

        CCStructure& cs = mesh->ccStructure("Tissue");

        for(auto c : cellAttr) {
            Tissue::CellData& cD = cellAttr[c.first];
            cD.resetMechanics(faceAttr);
        }
        for(CCIndex f : cs.faces()) {
            Tissue::FaceData& fD = faceAttr[f];
            fD.resetMechanics();
        }
        for(CCIndex e : cs.edges()) {
            Tissue::EdgeData& eD = edgeAttr[e];
            eD.resetMechanics();
        }
        /*for(CCIndex e : cs.vertices()) {
            Tissue::VertexData& vD = vMAttr[e];
            vD.resetMechanics();
        }*/
    }

    QString TissueName;
    QString TissueDualName;

    Process* rootProcess = 0;
    CellTissue cellTissue;
    std::map<IntIntPair, double> wallAreas;
    std::map<IntIntPair, std::set<CCIndex>> wallEdges;
    std::map<IntIntPair, std::map<string, double>> intercellularChems;
};

// Functions needed to serialized the structs attributes above,
// as the simple writechar cannot save more that primitive types
bool inline readAttr(Tissue::CellData& m, const QByteArray& ba, size_t& pos) {
    return m.read(ba, pos);
}
bool inline writeAttr(Tissue::CellData& m, QByteArray& ba) {
    return m.write(ba);
}

bool inline readAttr(Tissue::FaceData& m, const QByteArray& ba, size_t& pos) {
    return m.read(ba, pos);
}
bool inline writeAttr( Tissue::FaceData& m, QByteArray& ba) {
    return m.write(ba);
}

bool inline readAttr(Tissue::EdgeData& m, const QByteArray& ba, size_t& pos) {
    return m.read(ba, pos);
}
bool inline writeAttr(Tissue::EdgeData& m, QByteArray& ba) {
    return m.write(ba);
}

bool inline readAttr(Tissue::VertexData& m, const QByteArray& ba, size_t& pos) {
    return m.read(ba, pos);
}

bool inline writeAttr( Tissue::VertexData& m, QByteArray& ba) {
    return m.write(ba);
}


#endif // TISSUE_HPP
