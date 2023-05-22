#include "Root.hpp"
#include "PBD.hpp"

using namespace std;
using namespace ROOT;

bool Mechanics::initialize(QWidget* parent) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("MechanicalGrowth::step No current mesh"));

    if(!getProcess(parm("PBD Engine"), PBDProcess))
        throw(QString("Mechanics::initialize Cannot make PBD Engine") +
              parm("Mechanical Solver Process"));
    if(!getProcess(parm("Tissue Process"), tissueProcess))
        throw(QString("Chemicals::initialize Cannot make Tissue Process:") +
              parm("Tissue Process"));

    Verbose = parm("Verbose") == "True";

    if(Verbose)
        mdxInfo << "Mechanics::initialize" << endl;

    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    indexAttr = &mesh->indexAttr();

    Dt = parm("Dt").toDouble();

    PBDProcess->initialize(parent);

    wallStress = parm("Wall stress").toDouble();
    wallStressK1 = parm("Wall stress K1").toDouble();
    wallStressK2 = parm("Wall stress K2").toDouble();

    // External Forces
    QStringList list = parm("Gravity Direction").split(QRegExp(","));
    gravity[0] = list[0].toDouble();
    gravity[1] = list[1].toDouble();
    gravity[2] = list[2].toDouble();
    gravity *=  parm("Gravity Force").toDouble();

    // Hydrostatics
     for(uint i = 0; i < cellAttr.size(); i++) {
         auto it = cellAttr.begin();
         advance(it, i);
         Tissue::CellData& cD = it->second;
         if(cD.pressureMax == -1)
            cD.pressureMax = parm("Turgor Pressure").toDouble();
     }

    return true;
}

// basically it reset the chemicals
bool Mechanics::rewind(QWidget* parent) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("MechanicalGrowth::step No current mesh"));
    if(!mesh or mesh->file().isEmpty())
        throw(QString("Mechanics::rewind No current mesh, cannot rewind"));
    if(!tissueProcess)
        throw(QString("Mechanics::rewind tissue process not set"));
    if(Verbose)
        mdxInfo << "Mechanics::rewind" << endl;
    PBDProcess->initialize(parent);
    tissueProcess->resetMechanics();
    realTime = 0;
    userTime = 0;

    return false;
}

double Mechanics::calcNorm() {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("MechanicalGrowth::step No current mesh"));
    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");
    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    CCStructure& cs = mesh->ccStructure("Tissue");
    const CCIndexVec& vertices = cs.vertices();
    double avgNorm = 0.0, maxNorm = 0.0;
    //#pragma omp parallel for reduction(+ : avgNorm) reduction(max : maxNorm)
    for(uint i = 0; i < vertices.size(); i++) {
        CCIndex v = vertices[i];
        Tissue::VertexData& vD = vMAttr[v];
        CCIndexData& vI = indexAttr[v];
        double s = 0;
        s = norm(vI.pos - vD.prevPos) / Dt;
        avgNorm += s;
        if(maxNorm < s)
            maxNorm = s;
    }
    avgNorm /= vertices.size();
    double normal = 0.5 * (avgNorm + maxNorm);
    return normal;
}

bool Mechanics::step() {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("MechanicalGrowth::step No current mesh"));

    CCStructure& cs = mesh->ccStructure(mesh->ccName());
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    Tissue::FaceDataAttr& faceAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");


    // Internal Growing Forces on faces (replaced by PBD)
    for(uint i = 0; i < cs.faces().size(); i++) {
        CCIndex f = cs.faces()[i];
        Tissue::FaceData& fD = faceAttr[f];
        if (fD.type != Tissue::Membrane)
            continue;
        Tissue::CellData& cD = cellAttr[(*indexAttr)[f].label];

        cD.wallStress = wallStress;

        double auxinByArea = cD.auxin / cD.area;

        fD.stress = cD.wallStress *
                        ((pow(auxinByArea, 4)) / (pow(auxinByArea, 4) + pow(wallStressK1, 4))) *
                        ((pow(wallStressK2, 8)) / (pow(auxinByArea, 8) + pow(wallStressK2, 8)));

        fD.sigmaA = fD.stress *
                    (
                        (norm(fD.a1) * OuterProduct(fD.a1/norm(fD.a1), fD.a1/norm(fD.a1))) +
                        (norm(fD.a2) * OuterProduct(fD.a2/norm(fD.a2), fD.a2/norm(fD.a2)))
                    );

    }

    // Auxin relaxating effect on cell walls
    double baseWallEK = parm("Wall EK").toDouble();
    double auxinWallK1 = parm("Auxin-induced wall relaxation K1").toDouble();
    double auxinWallK2 = parm("Auxin-induced wall relaxation K2").toDouble();
    double pressureK = parm("Turgor Pressure Rate").toDouble();
    for(auto c : cellAttr) {
            Tissue::CellData& cD = cellAttr[c.first];
            double auxinByArea =  cD.auxin / cD.area;
            //cD.pressure = 1 + parm("Turgor Pressure").toDouble() * norm(cD.a1) * pow(auxinByArea, 4) / ( pow(wallStressK1, 4) +  pow(auxinByArea, 4)) ;
            cD.pressure += pressureK * Dt;
            if(cD.pressure > cD.pressureMax)
                cD.pressure = cD.pressureMax;
            for(CCIndex e : cD.perimeterEdges) {
                Tissue::EdgeData& eD = edgeAttr[e];
                if(eD.type == Tissue::Wall)
                  if(auxinWallK1 > 0)
                        eD.eStiffness = baseWallEK *  ( pow(auxinWallK1, 4) / ( pow(auxinWallK1, 4) +  pow(auxinByArea, 4)) +
                                                        pow(auxinByArea, 8) / ( pow(auxinWallK2, 8) +  pow(auxinByArea, 8))
                                                        );
            }

    }

    // Apply external/internal forces on vertices
    for(uint i = 0; i < cs.vertices().size(); i++) {
        CCIndex v = cs.vertices()[i];
        Tissue::VertexData& vD = vMAttr[v];
        calcForces(v, cs, *indexAttr,cellAttr, faceAttr, edgeAttr, vD);
    }

    PBDProcess->update(Dt);

    // the normal is our way to check for convergence
    prev_normal = normal;
    normal = calcNorm();

    // set user and real time
    userTime += Dt;
    if(normal - prev_normal > 0)
        realTime += Dt;

    // Root growth rate
    prevQCcentroid = QCcentroid;
    for(auto c : cellAttr)
        if(cellAttr[c.first].type == Tissue::QC)
            QCcentroid = cellAttr[c.first].centroid;
    growthRatesVector.push_back((QCcentroid - prevQCcentroid).norm() / Dt);
    if(growthRatesVector.size() > 30)
        growthRatesVector.erase(growthRatesVector.begin());

    if(++debug_step > parm("Debug Steps").toInt()) {
        double rootGR = 0;
        for(double gr : growthRatesVector)
            rootGR += gr;
        rootGR /= growthRatesVector.size();
        if(Verbose)
            mdxInfo << "Mechanics: " << "Norm: " << normal << " " << " realTime: " << realTime << " userTime: " << userTime
                << " growth rate " << rootGR <<  endl;
        debug_step = 0;

    }

    // check for convergence
    if(normal <= convergeThresh && convergenceLag > parm("Convergence Lag").toInt()) {
        convergenceLag = 0;
        return true;
    } else {
        convergenceLag += Dt;
        return false;
    }

}

Point3d Mechanics::calcForcesFace(CCIndex f,
    Tissue::CellData& cD, Tissue::FaceData& fD, Tissue::EdgeData& eD, Tissue::VertexData& vD ) {
    Point3d dx;

    // Hydrostatics
    /*
    if(eD.type == Tissue::Wall) {
        eD.pressureForce = cD.pressure * eD.outwardNormal[f] * eD.length;
        vD.forces.push_back(make_tuple(cD.label, "pressure", 0.5 * eD.pressureForce));
        dx += 0.5 * eD.pressureForce;
    }
    else
        eD.pressureForce = 0;
    */

    // Viscous Forces
    eD.sigmaForce = Point3d(0, 0, 0);
    // viscosity only active on wall edges
    if(eD.type == Tissue::Wall) {
            eD.sigmaForce = fD.sigmaA * eD.outwardNormal[f] * eD.length;
            vD.forces.push_back(make_tuple(cD.label, "sigmaAY", 0.5 * eD.sigmaForce));
            dx += 0.5 * eD.sigmaForce;
    }
    return dx;
}

// this function computes the force acting on the vtx v as member of the edge e.
// Its connected vtx n, with an eventual label to identify the cell exterting
// the force (0 if wall edge)
Point3d Mechanics::calcForcesEdge(
    const CCStructure& cs, const CCIndexDataAttr& indexAttr, CCIndex e, CCIndex v, int label) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("MechanicalGrowth::step No current mesh"));

    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");

    Tissue::EdgeData& eD = edgeAttr[e];
    Tissue::VertexData& vD = vMAttr[v];
    Point3d dx;
    auto eb = cs.edgeBounds(e);
    Point3d vPos = indexAttr[v].pos, nPos;
    if(eb.first == v)
        nPos = indexAttr[eb.second].pos;
    else
        nPos = indexAttr[eb.first].pos;

    // Mass Springs
    Point3d sigmaEeforce = eD.sigmaEe * ((vPos - nPos) / eD.length);
    vD.forces.push_back(make_tuple(label, "sigmaEe", sigmaEeforce));
    dx += sigmaEeforce;

    // Viscous Forces
    Point3d sigmaEvforce = eD.sigmaEv * ((vPos - nPos) / eD.length);
    vD.forces.push_back(make_tuple(label, "sigmaEv", sigmaEvforce));
    dx += sigmaEvforce;

    return dx;
}

void Mechanics::calcForces(CCIndex v, const CCStructure& cs, const CCIndexDataAttr& indexAttr, Tissue::CellDataAttr &cellAttr,
                           Tissue::FaceDataAttr &faceAttr,
                           Tissue::EdgeDataAttr& edgeAttr, Tissue::VertexData& vD) {
    vD.force = 0;
    vD.forces.clear();
    // Viscous forces on faces (replaced by PBD turgor + strain)
    CCIndexUSet edgeSet;
    if(vD.type == Tissue::Border) {
        for(auto flip : cs.matchV(v, CCIndex::Q, CCIndex::Q, CCIndex::Q)) {
            CCIndex f = flip.interior;
            if(f.isPseudocell() || indexAttr[f].label < 1)
                continue;
            Tissue::FaceData &fD = faceAttr[f];
            Tissue::EdgeData &eD0 = edgeAttr[flip.facet[0]];
            Tissue::EdgeData &eD1 = edgeAttr[flip.facet[1]];
            edgeSet.insert(flip.facet[0]);
            edgeSet.insert(flip.facet[1]);
            vD.force += calcForcesFace(f, cellAttr[indexAttr[f].label], fD, eD0, vD);
            vD.force += calcForcesFace(f, cellAttr[indexAttr[f].label], fD, eD1, vD);

        }

    }

    /* Elastic and viscous forces on edges (replaced by PBD distance constrain)
    for(CCIndex e : edgeSet) {
        int label = 0;
        if(vD.type == Tissue::Inside)
            for(CCIndex f : cs.incidentCells(v, 2))
                label = (indexAttr)[f].label;
        else if(vD.type == Tissue::Border)
            ;
        else
            throw(QString("Unknown Vertex Type: " + vD.type));
       vD.force += calcForcesEdge(cs, indexAttr, e, v, label);
    }*/


    // External Forces
    vD.force += gravity;
    //vD.forces.push_back(make_tuple(0, "gravity", gravity));
    if(cs.onBorder(v)) {
        Point3d friction = vD.velocity * parm("Friction").toDouble() ;
        vD.force -= friction;
        vD.forces.push_back(make_tuple(0, "friction", friction));
    }

    return;
}


bool MechanicalGrowth::initialize(QWidget* parent) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("MechanicalGrowth::step No current mesh"));
    if(!getProcess(parm("Mechanics Process"), mechanicsProcess))
        throw(QString("MechanicalGrowth::initialize Cannot make Mechanics:") + parm("Mechanics Process"));

    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
    Tissue::FaceDataAttr &faceAttr =
      mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");
    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    CCStructure& cs = mesh->ccStructure("Tissue");

    Verbose = parm("Verbose") == "True";

    if(Verbose)
        mdxInfo << "MechanicalGrowth::initialize" << endl;

    // strain constrain initial and restLength update
    for(uint i = 0; i < cellAttr.size(); i++) {
        auto it = cellAttr.begin();
        advance(it, i);
        Tissue::CellData& cD = it->second;
        if(cD.restArea <= 0) {
            cD.restArea = 0;
            for(CCIndex f : *cD.cellFaces)
                cD.restArea += indexAttr[f].measure;
        }
        if(cD.mfRORate == -1)
            cD.mfRORate = parm("MF reorientation rate").toDouble();
    }
    // init rest area
    for(CCIndex f : cs.faces()) {
        Tissue::FaceData& fD = faceAttr[f];
        if(fD.restAreaFace <= 0)
            fD.restAreaFace = indexAttr[f].measure;

    }
    // init rest length
    for(CCIndex e : cs.edges()) {
        Tissue::EdgeData& eD = edgeAttr[e];
        if(eD.restLength <= 0)
            eD.restLength = eD.length;
    }
    // int rest position
    for(CCIndex v : cs.vertices()) {
        Tissue::VertexData& vD = vMAttr[v];
        if(norm(vD.restPos) == 0)
            vD.restPos = indexAttr[v].pos;
    }

    return true;
}

bool MechanicalGrowth::step(double Dt) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("MechanicalGrowth::step No current mesh"));
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
    Tissue::FaceDataAttr &faceAttr =
      mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    CCStructure& cs = mesh->ccStructure("Tissue");



    /////////
    // Polarity
    ////////
    for(uint i = 0; i < cellAttr.size(); i++) {
        auto it = cellAttr.begin();
        advance(it, i);
        Tissue::CellData& cD = it->second;        
        if(cD.mfRORate != 0)
            cD.mfDelete = parm("MF Delete After Division") == "True";
        if(norm(cD.a1) == 0)
            cD.a1[1] = EPS;
        if(norm(cD.a2) == 0)
            cD.a2[0] = EPS;
        Point3d a1_inc, a2_inc;
        double threshold = parm("MF reorientation strainrate").toDouble();
        // Membrane stretch method
        if(parm("Polarity Method") == "Membrane Stretch") {
            Point3d stretch ;
            if(cD.area > cD.restArea) {
                for(CCIndex f : cD.perimeterFaces) {
                    Tissue::FaceData& fD = faceAttr[f];
                    if(fD.area > fD.restAreaFace) {
                        for(CCIndex e : cs.incidentCells(f, 1)) {
                            Tissue::EdgeData& eD = edgeAttr[e];
                            if(eD.length < eD.restLength) continue;
                            if(eD.type == Tissue::Wall) {
                                auto eb = cs.edgeBounds(e);
                                Point3d versor =
                                    indexAttr[eb.second].pos - indexAttr[eb.first].pos;
                                versor /= norm(versor);
                                if(versor * cD.a1 < 0)
                                    versor *= -1;
                                if(eD.strain > threshold)
                                    stretch += versor * eD.strain;
                            }
                        }
                    }
                }
            }
            a1_inc = stretch;
            a2_inc = -stretch;
        // Principal strains method
        } else if(parm("Polarity Method") == "Principal Strains") {
            if(parm("Strain Tensor") == "Green Strain Tensor") {
                if(cD.gMax > threshold)
                    a1_inc = Rotate(Point3d(cD.gMax,0,0), cD.gAngle);
                if(cD.gMin > threshold)
                    a2_inc = Rotate(Point3d(cD.gMin,0,0), cD.gAngle+(M_PI/2));
            } else if (parm("Strain Tensor") == "Shape Strain Tensor") {
                if(cD.sMax > threshold)
                    a1_inc = Rotate(Point3d(cD.sMax,0,0), cD.sAngle);
                if(cD.sMin > threshold)
                    a2_inc = Rotate(Point3d(cD.sMin,0,0), cD.sAngle+(M_PI/2));
            } else
                 throw(QString("Wrong Tensor"));
        // Cell Axis method
        } else if(parm("Polarity Method") == "Cell Axis") {
            if(cD.mfRORate > 0) {
                if(findClosestLineToLine(cD.axisMax,
                                         cD.a1, cD.a2) == cD.a1) {
                    cD.a1 = cD.axisMax/norm(cD.axisMax) * norm(cD.a1);
                    cD.a2 = cD.axisMin/norm(cD.axisMin) * norm(cD.a2);
                } else {
                    cD.a2 = cD.axisMax/norm(cD.axisMax) * norm(cD.a2);
                    cD.a1 = cD.axisMin/norm(cD.axisMin) * norm(cD.a1);
                }
                for(CCIndex e : cD.perimeterEdges) {
                    Tissue::EdgeData& eD = edgeAttr[e];
                    if(eD.strainRate < threshold)
                        continue;
                    auto eb = cs.edgeBounds(e);
                    Point3d versor = indexAttr[eb.first].pos - indexAttr[eb.second].pos;
                    versor /= norm(versor);
                    a1_inc += cD.a1/norm(cD.a1) * abs(cD.a1/norm(cD.a1) * versor * eD.strainRate);
                    a2_inc += cD.a2/norm(cD.a2) * abs(cD.a2/norm(cD.a2) * versor * eD.strainRate);
                }
            }
        }
        else
            throw(QString("Wrong Polarity Method"));

        // MF degradation
        Point3d a1U = norm(cD.a1) > 0 ? cD.a1 / norm(cD.a1) : Point3d(0, 0, 0);
        Point3d a2U = norm(cD.a2) > 0 ? cD.a2 / norm(cD.a2) : Point3d(0, 0, 0);
        Point3d degr_a1 =
            (cD.mfRORate > 0) ? parm("MF Degradation").toDouble() * a1U * norm(cD.a1): Point3d(0, 0, 0);
        Point3d degr_a2 =
            (cD.mfRORate > 0) ? parm("MF Degradation").toDouble() * a2U * norm(cD.a2) : Point3d(0, 0, 0);

        // MF final update
        cD.a1 += (a1_inc * cD.mfRORate - degr_a1) * Dt;
        cD.a2 += (a2_inc * cD.mfRORate - degr_a2) * Dt;

        // MF cannot be longer than 1
        if(norm(cD.a1) > 1)
            cD.a1 = shortenLength(cD.a1, norm(cD.a1) - 1);
        if(norm(cD.a2) > 1)
            cD.a2 = shortenLength(cD.a2, norm(cD.a2) - 1);
        // a1 is always the longest MF
        if(norm(cD.a1) < norm(cD.a2)) {
            Point3d tmp = cD.a1;
            cD.a1 = cD.a2;
            cD.a2 = tmp;
        }

        if(norm(cD.a1) == 0)
            cD.a1[1] = EPS;
        if(norm(cD.a2) == 0)
            cD.a2[0] = EPS;

        for(CCIndex f : (*cD.cellFaces)) {
            Tissue::FaceData& fD = faceAttr[f];
            fD.a1 = cD.a1;
            fD.a2 = cD.a2;
        }

        // division vector depends on the polarity method
        if(parm("Polarity Method") == "Membrane Stretch" ) {
            cD.divVector = cD.a1;
            // if we are in periclinal division rotate the vector 90 degrees
            if(!cD.periclinalDivision)
                cD.divVector = cD.divVector ^ Point3d(0,0,1);
        } else if (parm("Polarity Method") == "Cell Axis" || parm("Polarity Method") == "Principal Strains" ) {
            cD.divVector = cD.a1;
            if(!cD.periclinalDivision)
                cD.divVector = cD.a2;
        }

    }

    // Find the highest LRC cell (used later for zonation, to see how fare the cells are from the QC)
    double lrc = -BIG_VAL;
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        if(cD.type == Tissue::LRC && cD.centroid.y() > lrc)
            lrc = cD.centroid.y();
    }

    // Zonation and Resting values update
    double growthRateThresh = parm("Strain Threshold for Growth").toDouble();
    double wallsMaxGrowthRate = parm("Walls Growth Rate").toDouble();
    double elongationZone = parm("Elongation Zone").toDouble();
    double differentatiationZone = parm("Differentiation Zone").toDouble();
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        cD.lifeTime += Dt;
        // Zonation
        if(lrc != -BIG_VAL && elongationZone > 0 && differentatiationZone > 0 && elongationZone < differentatiationZone) {
            if(cD.centroid.y() - lrc >  elongationZone)
                cD.divisionAllowed = false;
            if(cD.centroid.y() - lrc >  differentatiationZone)
                cD.pressureMax = 1;
        }
        // Growth rates, rest lengths....
        // Disable growth update if this variable is zero, for debugging mostly
        if(cD.mfRORate == 0)
            continue;  ///
        // Disable growth update if cell division is not allowed (Crisanto's root) and we have reached max cell size
        if(!cD.divisionAllowed && cD.area > cD.cellMaxArea)
            continue;
        if(cD.growthRate > growthRateThresh)
            //cD.restArea += cD.area * areaMaxGrowthRate * Dt;
            cD.restArea = cD.area;
        for(CCIndex f : *cD.cellFaces) {
            Tissue::FaceData& fD = faceAttr[f];
            if(cD.growthRate > growthRateThresh)
                //fD.restAreaFace += fD.area * areaMaxGrowthRate * Dt;
                fD.restAreaFace = fD.area ;
            for(CCIndex e : cs.incidentCells(f, 1)) {
                Tissue::EdgeData& eD = edgeAttr[e];
                if(eD.type == Tissue::Wall) {
                    if(eD.strain >= growthRateThresh && cD.growthRate > growthRateThresh) {
                        //eD.restLength = eD.length;
                        double strainDiff = eD.strain - growthRateThresh;
                        eD.restLength += wallsMaxGrowthRate * strainDiff * Dt;
                        if( eD.restLength > eD.length) {
                            mdxInfo << "WARNING: restLength is updated inmediately" << endl;
                            eD.restLength = eD.length;
                        }
                    }
                } else if (eD.type == Tissue::Shear)
                    //if(eD.length > eD.restLength) // this seems to have very little effect
                    eD.restLength = eD.length ;
            }
        }
    }

    // synchonize columella initials and latep initials
    bool synchronized = true;
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        if(cD.type == Tissue::ColumellaInitial  || (cD.type == Tissue::EpLrcInitial)) {
            cD.divAlg = 2;
            if(cD.area < cD.cellMaxArea)
                synchronized = false;
        }
    }
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        if(cD.type == Tissue::ColumellaInitial || (cD.type == Tissue::EpLrcInitial  && cD.periclinalDivision)) {
            if(synchronized)
                cD.divisionAllowed = true;
            else
                cD.divisionAllowed = true;  ///////// disabled
        }
    }

    return false;
}

bool Chemicals::initialize(QWidget* parent) {
    mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("Chemicals::step No current mesh"));

    if(!getProcess(parm("Tissue Process"), tissueProcess))
        throw(QString("Chemicals::initialize Cannot make Tissue Process:") +
              parm("Tissue Process"));

    indexAttr = &mesh->indexAttr();
    cellAttr = &mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    faceAttr = &mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    vMAttr = &mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");

    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");

    // Set the auxin source from the top
    double auxin_source = parm("Auxin Source").toDouble();
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        if(cD.type == Tissue::Source  && auxin_source > 0)
            cD.auxinProdRate = auxin_source;
        if(cD.pinProdRate == -1)
            cD.pinProdRate = parm("Pin1 Basal Production Rate").toDouble();
        if(cD.pinInducedRate == -1)
            cD.pinInducedRate = parm("Pin1 Max Auxin-induced Expression").toDouble();
        if(cD.aux1ProdRate == -1)
            cD.aux1ProdRate = parm("Aux1 Basal Production Rate").toDouble();
        if(cD.aux1InducedRate == -1)
            cD.aux1InducedRate = parm("Aux1 Max Auxin-induced Expression").toDouble();
        if(cD.aux1MaxEdge == -1)
            cD.aux1MaxEdge = parm("Aux1 Max Amount Edge").toDouble();
    }

    Dt = parm("Dt").toDouble();
    convergeThresh = parm("Converge Threshold").toDouble();

    return true;
}

// basically it reset the chemicals
bool Chemicals::rewind(QWidget* parent) {
    mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("Chemicals::rewind No current mesh, cannot rewind"));
    tissueProcess->resetChems();
    userTime = 0;
    return false;
}


std::set<CCIndex> getBoundEdges(const CCStructure& cs, CCIndex e) {
    auto eb = cs.edgeBounds(e);
    std::set<CCIndex> edges;
    auto a = cs.incidentCells(eb.first, 1);
    auto b = cs.incidentCells(eb.second, 1);
    edges.insert(a.begin(),  a.end());
    edges.insert(b.begin(),  b.end());
    edges.erase(e);
    return (edges);
}

void Chemicals::calcDerivsEdge(const CCStructure& cs,
                               const CCIndexDataAttr& indexAttr,
                               Tissue::EdgeDataAttr& edgeAttr,
                               CCIndex e, double Dt) {
    Tissue::EdgeData& eD = edgeAttr[e];
    if(eD.type == Tissue::Wall) {
        double auxinDiffusionRate = parm("Auxin Intercellular Diffusion").toDouble();
        double kauxinDecay = parm("Auxin Decay Rate").toDouble();
        double kauxinMaxDecay = parm("Auxin Max Decay Rate").toDouble();
        double kauxinMaxEdge = parm("Auxin Max Amount Edge").toDouble();
        std::set<CCIndex> edges = getBoundEdges(cs, e);
        for(CCIndex en : edges) {
            Tissue::EdgeData& eDn = edgeAttr[en];
            if(eDn.type == Tissue::Wall) {
                // Auxin diffusion between adjacient intercellular edges
                double diffusion = auxinDiffusionRate * (eD.intercellularAuxin - eDn.intercellularAuxin) / (eD.length + eDn.length) * Dt;
                eD.intercellularAuxin -= diffusion;
                eDn.intercellularAuxin += diffusion;
            }
        }
        // auxin decay can be fast if the maximum amount is reached
        double edge_decay = eD.intercellularAuxin * (kauxinDecay + (kauxinMaxDecay - kauxinDecay) *
                                 (pow(eD.intercellularAuxin/eD.length,10) / (pow(kauxinMaxEdge,10) + pow(eD.intercellularAuxin/eD.length,10)))) * Dt;
        eD.intercellularAuxin -= edge_decay;
    }
}

// the TissueDual cell verteces
void Chemicals::calcDerivsCell(const CCStructure& cs,
                                  const CCStructure& csDual,
                                  const CCIndexDataAttr& indexAttr,
                                  Tissue::CellDataAttr& cellAttr,
                                  Tissue::EdgeDataAttr& edgeAttr,
                                  int label, double Dt) {


    Tissue::CellData &cD = cellAttr[label];
    if(cD.dualVertex.isPseudocell())
        return;

    // AUXIN derivatives
    cD.auxinFluxes.clear();
    double permeability = parm("Auxin Cell Permeability").toDouble();
    double Kaux1 = parm("Aux1-auxin import rate").toDouble();
    double Kpin1 = parm("Pin1-auxin export rate").toDouble();
    double kauxinMaxCell = parm("Auxin Max Amount Cell").toDouble();
    if(cD.type == Tissue::Source)
        kauxinMaxCell = BIG_VAL;
    double kauxinMaxEdge = parm("Auxin Max Amount Edge").toDouble();
    double kauxinDecay = parm("Auxin Decay Rate").toDouble();
    double kauxinMaxDecay = parm("Auxin Max Decay Rate").toDouble();
    double diffusion = 0, activeTransport = 0;
    for(CCIndex vn : csDual.neighbors(cD.dualVertex)) {
        // changes relative to the current neighbour
        int labeln = indexAttr[vn].label;
        double cell_diffusion = 0, cell_activeTransport = 0;
        cD.auxinFluxes[labeln] = 0;
        for(CCIndex e : tissueProcess->wallEdges[std::make_pair(label, labeln)]) {
            Tissue::EdgeData& eD = edgeAttr[e];
            double edge_diffusion = permeability * (eD.intercellularAuxin - cD.auxin) * 1. / cD.area * eD.length;
            double edge_import = Kaux1 * eD.intercellularAuxin * eD.Aux1[label]  * 1. / 1.; // we assume intercell space 1 nm thick
            double edge_export = Kpin1 * cD.auxin * eD.Pin1[label] * 1. / cD.area; // FIXME, is this wrong? where is the length?
            double edge_transport = edge_import - edge_export;
            // check for maximums
            if(edge_diffusion > 0)
                edge_diffusion *= pow(kauxinMaxCell,10) / (pow(kauxinMaxCell,10) + pow(cD.auxin/cD.area,10));
            if(edge_transport > 0)
                edge_transport *= pow(kauxinMaxCell,10) / (pow(kauxinMaxCell,10) + pow(cD.auxin/cD.area,10));
            if(edge_diffusion < 0)
                edge_diffusion *= pow(kauxinMaxEdge,10) / (pow(kauxinMaxEdge,10) + pow(eD.intercellularAuxin/eD.length,10));
            if(edge_transport < 0)
                edge_transport *= pow(kauxinMaxEdge,10) / (pow(kauxinMaxEdge,10) + pow(eD.intercellularAuxin/eD.length,10));

            eD.intercellularAuxin += (- edge_diffusion - edge_transport) * Dt;
            eD.auxinRatio[label] = edge_transport /* 1. / eD.length*/ * Dt;
            cD.auxinFluxes[labeln] += (cD.centroid - eD.midPoint)/norm(cD.centroid - eD.midPoint) * edge_transport * Dt;
            cell_diffusion += edge_diffusion;
            cell_activeTransport += edge_transport;
        }
        // update the total changes of the cell
        diffusion += cell_diffusion;
        activeTransport += cell_activeTransport;
    }
    double Kbase = parm("Auxin Basal Production Rate").toDouble();
    double basalProduction = Kbase + cD.auxinProdRate * (pow(kauxinMaxCell,10) / (pow(kauxinMaxCell,10) + pow(cD.auxin/cD.area,10)));
    double decay = cD.auxin * (kauxinDecay + (kauxinMaxDecay - kauxinDecay) *
                              (pow(cD.auxin/cD.area,10) / (pow(kauxinMaxCell,10) + pow(cD.auxin/cD.area,10))));
    cD.prevAuxin = cD.auxin;
    cD.auxin += (basalProduction - decay + diffusion + activeTransport) * Dt;
    if(cD.selected)
        mdxInfo << "basalProduction: " << basalProduction << " decay: " << decay << " diffusion " << diffusion << " activeTransport: " << activeTransport << endl;
    debugs["Average Auxin Diffusion"] += -diffusion * Dt;
    debugs["Average Auxin Export"] += activeTransport * Dt;

    // AUX1 Cytoplasmic
    double AUX1basalProduction = cD.aux1ProdRate;
    double AUX1reg = cD.aux1InducedRate;
    double AUX1Kaux = parm("Aux1 Half-max Auxin-induced K").toDouble();
    double AUX1max = parm("Aux1 Max Concentration").toDouble();
    double AUX1decayRate = parm("Aux1 Cytoplasmic Decay Rate").toDouble();
    double AUX1inducedExpression = AUX1reg *
                                    (
                                     (pow(cD.auxin / cD.area, 2) / (pow(AUX1Kaux, 2) + pow(cD.auxin / cD.area, 2))) *
                                     (pow(AUX1max, 8) / (pow(AUX1max, 8) + pow(cD.Aux1 / cD.area, 8)))
                                    );
    double AUX1decay = AUX1decayRate * cD.Aux1;
    cD.Aux1 += (AUX1basalProduction + AUX1inducedExpression  - AUX1decay) * Dt;
    if(cD.type == Tissue::Source)
        cD.Aux1 = 0;
    if(cD.type == Tissue::Substrate) // Substrate has no auxin so is does not express AUX1 by itself
        cD.Aux1 = 100;
    debugs["Average AUX1 Expression"] += AUX1inducedExpression * Dt;

    // AUX1 membranes derivatives
    double AUX1Krate = parm("Aux1 Max Trafficking Rate").toDouble();
    double AUX1MaxEdge = cD.aux1MaxEdge;
    for(CCIndex vn : csDual.neighbors(cD.dualVertex)) {
        for(CCIndex e : tissueProcess->wallEdges[std::make_pair(label, indexAttr[vn].label)]) {
    //for(CCIndex e : cD.perimeterEdges){{
            Tissue::EdgeData& eD = edgeAttr[e];
            double Aux1Sensitivity = eD.length / cD.perimeter;
            double traffickedAUX1 = AUX1Krate * Aux1Sensitivity * cD.Aux1
                     * (pow(AUX1MaxEdge,10) / (pow(AUX1MaxEdge,10) + pow(eD.Aux1[label]/eD.length,10)));
            double decayedAUX1 = eD.Aux1[label] * AUX1decayRate;
            eD.Aux1[label] += (traffickedAUX1 - decayedAUX1) * Dt;
            cD.Aux1 -= traffickedAUX1 * Dt;
            debugs["Average AUX1 Trafficking"] += traffickedAUX1 * Dt;
        }
    }


    // Pin1 cytoplasmic derivatives
    double pinBasalProduction = cD.pinProdRate;
    double Preg = cD.pinInducedRate;
    double Kaux = parm("Pin1 Half-max Auxin-induced K").toDouble();
    double Pmax = parm("Pin1 Max Concentration").toDouble();
    double inducedExpression = Preg *
                                    (
                                     (pow(cD.auxin / cD.area, 2) / (pow(Kaux, 2) + pow(cD.auxin / cD.area, 2))) *
                                     (pow(Pmax, 8) / (pow(Pmax, 8) + pow(cD.Pin1 / cD.area, 8)))
                                    );
    /* High auxin induces PIN decay
    double dauxin = parm("Pin1 Max Auxin-induced Decay").toDouble();
    double Kdegr = parm("Pin1 Half-max Auxin-induced Decay").toDouble();
    double inducedDecay = dauxin * (pow(cD.auxin, 2) / (pow(Kdegr, 2) + pow(cD.auxin, 2))) * cD.Pin1;
    */
    double inducedDecay = 0;
    decay = parm("Pin1 Cytoplasmic Decay Rate").toDouble() * cD.Pin1;
    cD.Pin1 += (pinBasalProduction + inducedExpression - inducedDecay - decay) * Dt;
    debugs["Average PIN1 Expression"] += inducedExpression * Dt;

    // Pin1 membranes derivatives
    double Krate = parm("Pin1 Max Trafficking Rate").toDouble() * (cD.auxin / (cD.auxin + 1)); //// FIXME no coefficient here
    double pinCytDecayRate = parm("Pin1 Cytoplasmic Decay Rate").toDouble();
    double pinMemDecayRate = parm("Pin1 Membrane Max Decay Rate").toDouble();
    double pinMaxEdge = parm("Pin1 Max Amount Edge").toDouble();
    for(CCIndex vn : csDual.neighbors(cD.dualVertex)) {
        for(CCIndex e : tissueProcess->wallEdges[std::make_pair(label, indexAttr[vn].label)]) {
    //for(CCIndex e : cD.perimeterEdges){{
            Tissue::EdgeData& eD = edgeAttr[e];
            double traffickedPin1 = Krate * eD.pin1Sensitivity[label] * cD.Pin1
                     * (pow(pinMaxEdge,10) / (pow(pinMaxEdge,10) + pow(eD.Pin1[label]/eD.length,10)));
            // higher degradation on edges with low sensitivity or low cell auxin, to allow faster turnout
            double decayedPin1 = eD.Pin1[label]
                    *  (pinCytDecayRate + (pinMemDecayRate - pinCytDecayRate) *
                        (                            
                            (0.1 + 1)   //// this should be fixed
                          + (pow(0.01, 10) / (pow(eD.pin1Sensitivity[label], 10) +  pow(0.01, 10))) //// FIXME no regulated coefficient 0.01 here
                        )
                       );
            eD.Pin1[label] += (traffickedPin1 - decayedPin1) * Dt;
            cD.Pin1 -= traffickedPin1 * Dt;
            debugs["Average PIN1 Trafficking"] += traffickedPin1 * Dt;
        }
    }

    // PHOSPHO cytoplasmic derivatives
    double PINOIDbasalProductionRate = parm("PINOID Basal Production Rate").toDouble();
    double PINOIDdecayRate = parm("PINOID Decay Rate").toDouble();
    double PP2AbasalProductionRate = parm("PP2A Basal Production Rate").toDouble() ;
    double PP2AdecayRate  = parm("PP2A Decay Rate").toDouble();

    // PHOSPHO membranes derivatives
    double PINOIDKrate = parm("PINOID Trafficking Rate").toDouble();
    double PINOIDMaxEdge = parm("PINOID Max Amount Edge").toDouble();
    double PP2AKrate = parm("PP2A Trafficking Rate").toDouble();
    double PP2AMaxEdge = parm("PP2A Max Amount Edge").toDouble();
    double PINOID_disp_K = parm("PINOID Displacing K").toDouble();
    double PINOID_fluidity_K = parm("PINOID Fluidity K").toDouble();
    double AUXIN_PP2A_Relief_T = parm("Auxin-PP2A Relief T").toDouble();
    double AUXIN_PP2A_Relief_K = parm("Auxin-PP2A Relief K").toDouble();
    double PINOID_PP2A_Trafficking_Toggle_K = parm("PINOID-PP2A Trafficking Toggle K").toDouble();
    double PINOID_PP2A_Disassociation_Toggle_K = parm("PINOID-PP2A Disassociation Toggle K").toDouble();
    double Geom_PP2A_Relief_T = parm("Geom-PP2A Relief T").toDouble();
    double Geom_PP2A_Relief_K = parm("Geom-PP2A Relief K").toDouble();
    double MF_PP2A_Relief_T = parm("MF-PP2A Relief T").toDouble();
    double MF_PP2A_Relief_K = parm("MF-PP2A Relief K").toDouble();
    double PINOIDdilutionRate = parm("PINOID Dilution Rate").toDouble();
    double PP2AdilutionRate = parm("PP2A Dilution Rate").toDouble();

    double trafficked_PINOID_total = 0, trafficked_PP2A_total = 0;
    double disassociated_PINOID_total = 0, disassociated_PP2A_total = 0;
    /*for(CCIndex vn : csDual.neighbors(cD.dualVertex)) {
        for(CCIndex e : tissueProcess->wallEdges[std::make_pair(label, indexAttr[vn].label)]) {*/
    for(CCIndex e : cD.perimeterEdges){{
            Tissue::EdgeData& eD = edgeAttr[e];
            double PINOID_conc_eD = eD.PINOID[label] / eD.length;
            double PP2A_conc_eD = eD.PP2A[label] / eD.length;
            double auxin_ratio = eD.auxinGrad[label] > 0 ? eD.auxinGrad[label] : 0;
            // exocytosis
            double trafficked_PINOID = Dt * cD.PINOID
                        * PINOIDKrate * (eD.length / cD.perimeter) // equal trafficking for edge
                        * (pow(cD.auxin / cD.area, 2) / (pow(0.01, 2) + pow(cD.auxin / cD.area, 2))) // auxin promotes trafficking
                        * (pow(PINOIDMaxEdge,10) / (pow(PINOIDMaxEdge,10) + pow(PINOID_conc_eD,10))) // Max amount on edge
                        ;
            double trafficked_PP2A = Dt * cD.PP2A * (eD.length / cD.perimeter)
                        * (  // not equal trafficking for edge
                           PP2AKrate +
                           AUXIN_PP2A_Relief_T * (pow(auxin_ratio, 4) / (pow(auxin_ratio, 4) + pow(AUXIN_PP2A_Relief_K, 4))) +
                           Geom_PP2A_Relief_T * (pow(1-eD.geomImpact[label], 4) / (pow(1-eD.geomImpact[label], 4) + pow(1-Geom_PP2A_Relief_K, 4))) +
                           MF_PP2A_Relief_T * (pow(1-eD.MFImpact[label], 4) / (pow(1-eD.MFImpact[label], 4) + pow(1-MF_PP2A_Relief_K, 4)))
                          )
                        * (pow(PINOID_PP2A_Trafficking_Toggle_K, 8) / (pow(PINOID_PP2A_Trafficking_Toggle_K, 8) + pow(PINOID_conc_eD, 8))) // not used
                        * (pow(cD.auxin / cD.area, 2) / (pow(0.01, 2) + pow(cD.auxin / cD.area, 2))) // auxin promotes trafficking
                        * (pow(PP2AMaxEdge,10) / (pow(PP2AMaxEdge,10) + pow(PP2A_conc_eD,10))) // Max amount on edge
                      ;
            trafficked_PINOID_total += trafficked_PINOID;
            trafficked_PP2A_total += trafficked_PP2A;

            double disassociated_PINOID = 0;
            double disassociated_PP2A = Dt * eD.PP2A[label]
                        * (pow(PINOID_conc_eD, 8) / (pow(PINOID_PP2A_Disassociation_Toggle_K, 8) + pow(PINOID_conc_eD, 8)))
                      ;
            disassociated_PINOID_total += disassociated_PINOID;
            disassociated_PP2A_total += disassociated_PP2A;

            std::vector<CCIndex>::iterator it = std::find(cD.perimeterEdges.begin(), cD.perimeterEdges.end(), e);
            int curr_index = std::distance(cD.perimeterEdges.begin(), it);
            int next_index = curr_index + 1;
            int prev_index = curr_index - 1;
            if((uint)next_index == cD.perimeterEdges.size())
                next_index = 0;
            if(prev_index == -1)
                prev_index = cD.perimeterEdges.size() - 1;
            CCIndex en = cD.perimeterEdges[next_index];
            CCIndex ep = cD.perimeterEdges[prev_index];
            Tissue::EdgeData& eDn = edgeAttr[en];
            Tissue::EdgeData& eDp = edgeAttr[ep];
            double PP2A_conc_eDn = eDn.PP2A[label] / eDn.length;
            double PP2A_conc_eDp = eDp.PP2A[label] / eDp.length;
            double PINOID_conc_eDn = eDn.PINOID[label] / eDn.length;
            double PINOID_conc_eDp = eDp.PINOID[label] / eDp.length;

            // Dilution
            double dilution_PINOID_n = 0;
            double dilution_PP2A_n = 0;
            double dilution_PINOID_p = 0;
            double dilution_PP2A_p = 0;
            //if(!cs.onBorder(en))
            {
                 dilution_PINOID_n = PINOIDdilutionRate * (PINOID_conc_eD - PINOID_conc_eDn)  * Dt;
                 dilution_PP2A_n = PP2AdilutionRate * (PP2A_conc_eD - PP2A_conc_eDn)  * Dt;
            }
            //if(!cs.onBorder(ep))
            {
                dilution_PINOID_p = PINOIDdilutionRate * (PINOID_conc_eD - PINOID_conc_eDp)  * Dt;
                dilution_PP2A_p = PP2AdilutionRate * (PP2A_conc_eD - PP2A_conc_eDp)  * Dt;
            }
            eD.PP2A[label] -= dilution_PP2A_n + dilution_PP2A_p;
            eD.PINOID[label] -=  dilution_PINOID_n + dilution_PINOID_p;
            eDn.PP2A[label] += dilution_PP2A_n;
            eDn.PINOID[label] += dilution_PINOID_n;
            eDp.PP2A[label] += dilution_PP2A_p;
            eDp.PINOID[label] += dilution_PINOID_p;

            CCIndex ec;
            /*if(cs.onBorder(en))
                ec = ep;
            else if (cs.onBorder(ep))
                ec = en;
            else*/
                ec = (rand() > RAND_MAX/2) ? en : ep;
            Tissue::EdgeData& eDc = edgeAttr[ec];

            double PP2A_conc_eDc = eDc.PP2A[label] / eDc.length;
            double PINOID_conc_eDc = eDc.PINOID[label] / eDc.length;

            // PINOIDS moving toward the lowest
            if(abs(PP2A_conc_eDc - PP2A_conc_eD ) > EPS){
                if(PP2A_conc_eDc <= PP2A_conc_eD || rand() < PINOID_fluidity_K * RAND_MAX * (PP2A_conc_eD / (PP2A_conc_eDc + EPS))) {
                      double disp_PINOID = PINOID_disp_K / (eD.length + eDc.length) * eD.PINOID[label] * Dt;
                      disp_PINOID *= (pow(PINOIDMaxEdge,10) / (pow(PINOIDMaxEdge,10) + pow(PINOID_conc_eDc,10))); // Max amount on edge
                      eDc.PINOID[label] += disp_PINOID;
                      eD.PINOID[label] -= disp_PINOID;
                }
            }

            // decay
            double decayedPINOID = Dt * eD.PINOID[label] * PINOIDdecayRate;
            double decayedPP2A = Dt * eD.PP2A[label] * PP2AdecayRate;
            eD.PP2A[label] += trafficked_PP2A - disassociated_PP2A - decayedPP2A;
            eD.PINOID[label] += trafficked_PINOID -  disassociated_PINOID - decayedPINOID;

        }
    }
    double PINOIDbasalProduction = PINOIDbasalProductionRate * Dt;
    double PP2AbasalProduction = PP2AbasalProductionRate * Dt;
    double PINOIDdecay = PINOIDdecayRate * cD.PINOID * Dt;
    double PP2Adecay = PP2AdecayRate * cD.PP2A * Dt;
    cD.PINOID += PINOIDbasalProduction - trafficked_PINOID_total + disassociated_PINOID_total - PINOIDdecay;
    cD.PP2A += PP2AbasalProduction - trafficked_PP2A_total + disassociated_PP2A_total - PP2Adecay;

    // Division Inhibitor and Promoter
    double divInhibitorPermeability = parm("Division Inhibitor Permeability").toDouble();
    double divPromoterPermeability = parm("Division Promoter Permeability").toDouble();
    for(CCIndex vn : csDual.neighbors(cD.dualVertex)) {
        int labeln = indexAttr[vn].label;
        Tissue::CellData& cDn = cellAttr[labeln];
        for(CCIndex e : tissueProcess->wallEdges[std::make_pair(label, labeln)]) {
            Tissue::EdgeData& eD = edgeAttr[e];
            double inhibitorDiffusion = 0.5 * divInhibitorPermeability * eD.length * (cDn.divInhibitor/cDn.area - cD.divInhibitor/cD.area)  ;
            cD.divInhibitor += inhibitorDiffusion * Dt;
            cDn.divInhibitor -= inhibitorDiffusion * Dt;
            double promoterDiffusion = 0.5 * divPromoterPermeability * eD.length * (cDn.divPromoter/cDn.area - cD.divPromoter/cD.area)   ;
            cD.divPromoter += promoterDiffusion * Dt;
            cDn.divPromoter -= promoterDiffusion * Dt;
        }
    }
    // Promoter
    double divBase = parm("Division Promoter Basal Production Rate").toDouble();
    double divKmax = parm("Division Promoter Max Auxin-induced Expression").toDouble();
    double divK = parm("Division Promoter Half-max Auxin-induced K").toDouble();
    int divN = parm("Division Promoter Half-max Auxin-induced n").toInt();
    double divInduced = 0;
    if(cD.divPromoter/cD.area < 5)
        if(cD.type == Tissue::QC || cD.type == Tissue::ColumellaInitial || cD.type == Tissue::VascularInitial || cD.type == Tissue::CEI || cD.type == Tissue::CEID || cD.type == Tissue::Columella || cD.type == Tissue::EpLrcInitial)
            divInduced = divKmax *
                                    (
                                     (pow(cD.auxin / cD.area, divN) / (pow(divK, divN) + pow(cD.auxin / cD.area, divN)))
                                     //(pow(AUX1max, 8) / (pow(AUX1max, 8) + pow(cD.Aux1 / cD.area, 8)))
                                    );
    double divDecay = cD.divPromoter * parm("Division Promoter Decay Rate").toDouble();
    cD.divPromoter += (divBase - divDecay + divInduced) * Dt;
    // Inhibitor
    divBase = parm("Division Inhibitor Basal Production Rate").toDouble();
    divKmax = parm("Division Inhibitor Max Promoter-induced Expression").toDouble();
    divK = parm("Division Inhibitor Half-max Promoter-induced K").toDouble();
    divN = parm("Division Inhibitor Half-max Promoter-induced n").toInt();
    divInduced = 0;
    if(cD.divInhibitor/cD.area < 5)
        //if(cD.type == Tissue::QC || cD.type == Tissue::ColumellaInitial || cD.type == Tissue::VascularInitial || cD.type == Tissue::CEI || cD.type == Tissue::CEID || cD.type == Tissue::Columella || cD.type == Tissue::EpLrcInitial)
            divInduced = divKmax *
                                    (
                                     (pow(cD.divPromoter / cD.area, divN) / (pow(divK, divN) + pow(cD.divPromoter / cD.area, divN)))
                                     //(pow(AUX1max, 8) / (pow(AUX1max, 8) + pow(cD.Aux1 / cD.area, 8)))
                                    );
    divDecay = cD.divInhibitor * parm("Division Inhibitor Decay Rate").toDouble();
    cD.divInhibitor += (divBase - divDecay + divInduced) * Dt;



    if(cD.selected) {
        for(CCIndex f : *cD.cellFaces)
          if(indexAttr[f].selected)
              for(CCIndex ez : cs.incidentCells(f, 1)){
                  Tissue::EdgeData& eD = edgeAttr[ez];
                    cout << cD.label << endl;
                    cout <<  "User Time: " << userTime <<  " PIN mem " << ez <<  " " << eD.Pin1[label] << " AUX1 mem " << ez <<  " " << eD.Aux1[label]  <<  endl;

                }
    }

    // Undefined are like dead cells
    if(cD.type == Tissue::Undefined) {
        cD.auxin = 0;
        cD.Aux1 = 0;
        cD.Pin1 = 0;
        cD.PINOID = 0;
        cD.PP2A = 0;
        cD.auxinProdRate = 0;
        cD.pinProdRate = 0;
        cD.pinInducedRate = 0;
    }

}


double Chemicals::calcNorm() {
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    double avgNorm = 0.0, maxNorm = 0.0;
    //#pragma omp parallel for reduction(+ : avgNorm) reduction(max : maxNorm)
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        double s = 0;
        s = norm(cD.auxin - cD.prevAuxin) / Dt;
        avgNorm += s;
        if(maxNorm < s)
            maxNorm = s;
    }
    avgNorm /= cellAttr.size();
    double normal = 0.5 * (avgNorm + maxNorm);
    return normal;
}

bool Chemicals::update() {

    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("Chemicals::step No current mesh"));
    if(!tissueProcess)
        throw(QString("Chemicals::step No Tissue Process"));

    CCIndexDataAttr indexAttr = mesh->indexAttr();
    CCStructure& cs = mesh->ccStructure("Tissue");
    CCStructure& csDual = mesh->ccStructure("TissueDual");
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");

    // some debug infos
    debugs.clear();

    // coordinates of the root tip
    Point3d root_tip;
    root_tip = BIG_VAL;
    for(uint i = 0; i < cellAttr.size(); i++) {
        auto it = cellAttr.begin();
        advance(it, i);
        Tissue::CellData& cD = it->second;
        if(cD.centroid.y() < root_tip.y())
            root_tip = cD.centroid;
    }

    // Calculate chemicals changes on cells
    for(auto c : cellAttr)
            calcDerivsCell(cs, csDual, indexAttr, cellAttr, edgeAttr, c.first, Dt);

    // Calculate chemicals changes on edges
    for(CCIndex e : cs.edges())
            calcDerivsEdge(cs, indexAttr, edgeAttr, e, Dt);

    // compute the pseudo-gradient of diffusion, the auxin flux vector and PIN sensitivity
    //#pragma omp parallel for
    for(uint i = 0; i < cellAttr.size(); i++) {
        auto it = cellAttr.begin();
        advance(it, i);
        Tissue::CellData& cD = it->second;
        int label = cD.label;
        if(cD.dualVertex.isPseudocell())
            continue;
        // calculate auxin gradient between adjacient cells
        cD.auxinFluxVector = Point3d(0., 0., 0.);
        for(CCIndex vn : csDual.neighbors(cD.dualVertex)) {
            int labeln = indexAttr[vn].label;
            if(cellAttr.find(labeln) ==  cellAttr.end())
                continue;
            cD.auxinFluxVector += cD.auxinFluxes[labeln];
        }
         // Auxin ratio gradients (used by the polarizer model)
        for(CCIndex e : cD.perimeterEdges) {
            Tissue::EdgeData& eD = edgeAttr[e];
            eD.auxinGrad[label] = 0;
            eD.auxinGrad[label] = eD.auxinRatio[label];
            std::vector<CCIndex>::iterator it = std::find(cD.perimeterEdges.begin(), cD.perimeterEdges.end(), e);
            int curr_index = std::distance(cD.perimeterEdges.begin(), it);
            double distance = 0;
            for(uint index = curr_index + 1; index <= cD.perimeterEdges.size() - 1; index++) {
                Tissue::EdgeData& eDn = edgeAttr[cD.perimeterEdges[index]];
                distance += norm(eDn.midPoint - eD.midPoint);
                eD.auxinGrad[label] +=  eDn.auxinRatio[label]  * 1. / distance;
            }
            distance = 0;
            for(int index = curr_index - 1; index >= 0 ; index--) {
                Tissue::EdgeData& eDn = edgeAttr[cD.perimeterEdges[index]];
                distance += norm(eDn.midPoint - eD.midPoint);
                 eD.auxinGrad[label] +=  eDn.auxinRatio[label] * 1. / distance;
            }
            eD.auxinGrad[label] *= 1 / cD.area * 1000;
        }
        // calculate PIN sensitivity
        double kauxinMaxCell = parm("Auxin Max Amount Cell").toDouble();
        double KauxinFlux = parm("Auxin-Flux Impact Half-max").toDouble();
        double Kpinoid = parm("PINOID Impact Half-max").toDouble();
        double Kmf = parm("MF Impact Half-max").toDouble();
        double Kgeom = parm("Geometry Impact Half-max").toDouble();
        double KpinSuppAuxin = parm("Pin1 Sensitivity Suppression by Auxin Amount").toDouble();
        double KPin1MF = parm("Pin1 Sensitivity MF K").toDouble();
        double KPin1auxinflux = parm("Pin1 Sensitivity Auxin-flux K").toDouble();
        double KPin1geom = parm("Pin1 Sensitivity Geometry K").toDouble();
        double KPin1MFauxinflux = parm("Pin1 Sensitivity MF+Auxin-flux K").toDouble();
        double KPin1MFgeom  = parm("Pin1 Sensitivity MF+Geometry K").toDouble();
        double KPin1geomauxinflux  = parm("Pin1 Sensitivity Auxin-flux+Geometry K").toDouble();
        double KPINOIDMax = parm("PINOID Max Amount Edge").toDouble();
        double normSum = 0;
        for(CCIndex e : cD.perimeterEdges) {
            Tissue::EdgeData& eD = edgeAttr[e];
            for(CCIndex f : cs.incidentCells(e, 2))
                if(indexAttr[f].label == label) {
                    // Calculate auxin flux impact
                    if(parm("Auxin Polarity Method") == "Flow") {
                        eD.auxinFluxImpact[label] = 0;
                        double auxinImpactValue = 0;
                        if(norm(cD.auxinFluxVector) > 0)
                            auxinImpactValue =
                                norm(cD.auxinFluxVector) * cD.auxinFluxVector/norm(cD.auxinFluxVector) * eD.outwardNormal[f] * eD.length;
                        if(auxinImpactValue < 0 || cD.auxin/cD.area > kauxinMaxCell)
                                if(parm("Pin1 Sensitivity Suppression by Auxin Max Cell") == "True")
                                    auxinImpactValue = 0;
                        eD.auxinFluxImpact[label] =
                            (pow(auxinImpactValue, 4) / (pow(KauxinFlux, 4) + pow(auxinImpactValue, 4))) ;
                    }
                    else if (parm("Auxin Polarity Method") == "PINOID")
                         eD.auxinFluxImpact[label] =
                             (pow(eD.PINOID[label]/(KPINOIDMax*eD.length), 4) / (pow(Kpinoid, 4) + pow(eD.PINOID[label]/(KPINOIDMax*eD.length), 4))) ;
                    else
                        throw(QString("Choose a correct auxin flux method, you ass"));
                    // Calculate MF impact
                    double MFvalue = cD.a1/norm(cD.a1) * eD.outwardNormal[f] ;
                    if(MFvalue < 0)
                        MFvalue *= -1;
                    eD.MFImpact[label] =  norm(cD.a1) * // eD.length *
                        (pow(MFvalue, 4) / (pow(Kmf, 4) + pow(MFvalue, 4))) ;
                    if(cs.onBorder(e))
                        eD.MFImpact[label] = 0;
                    // Geometry
                    double geomRatio = 1. - norm(cD.axisMin) / norm(cD.axisMax);
                    double geomAlign = cD.axisMax/norm(cD.axisMax) * eD.outwardNormal[f];
                    eD.geomImpact[label] = geomAlign * (pow(geomRatio, 4) / (pow(Kgeom, 4) + pow(geomRatio, 4))) ;
                    if(eD.geomImpact[label] < 0)
                        eD.geomImpact[label] *= -1;
                    // Calculate raw Pin sensitivity
                    eD.pin1SensitivityRaw[label] =
                              //eD.length *
                              (pow(KpinSuppAuxin, 10) / (pow(cD.auxin/cD.area, 10) + pow(KpinSuppAuxin, 10)))
                            *   (
                                (KPin1MF * eD.MFImpact[label])
                            +   (KPin1auxinflux * eD.auxinFluxImpact[label])
                            +   (KPin1geom * eD.geomImpact[label])
                            +   (KPin1MFauxinflux * eD.MFImpact[label] * eD.auxinFluxImpact[label])
                            +   (KPin1MFgeom * eD.MFImpact[label] * eD.geomImpact[label])
                            +   (KPin1geomauxinflux * eD.auxinFluxImpact[label] * eD.geomImpact[label])
                               );
                            ;
                    // Columella has PIN3, uniform pins all around
                    if(cs.onBorder(e) || (cD.type == Tissue::Columella && parm("Columella Auto-Efflux") == "True"))
                       eD.pin1SensitivityRaw[label] = 0;
                    // simulate PIN4, does not work if the root bends
                    if(parm("Simulate PIN4") == "True")
                        for(CCIndex vn : csDual.neighbors(cD.dualVertex)) {
                            int labeln = indexAttr[vn].label;
                            if(cellAttr.find(labeln) ==  cellAttr.end())
                                continue;
                            Tissue::CellData& cDn = cellAttr[labeln];
                            if(cDn.type == Tissue::QC)
                                if(norm(eD.midPoint - root_tip) > norm(cD.centroid - root_tip))
                                    eD.pin1SensitivityRaw[label] = 0;
                        }
                    // Sources have stronger directional flow
                    if(cD.type == Tissue::Source)
                       eD.pin1SensitivityRaw[label] = 10 * eD.MFImpact[label];
                    if(parm("Pin1 Sensitivity Average Method") == "Soft-max")
                        normSum += exp(eD.pin1SensitivityRaw[label]);
                    else if(parm("Pin1 Sensitivity Average Method") == "Arithmetic Average")
                        normSum += eD.pin1SensitivityRaw[label];
                }
        }
        // now softmax the Pin sensitivity on the edges
        for(CCIndex e : cD.perimeterEdges) {
            Tissue::EdgeData& eD = edgeAttr[e];
            for(CCIndex f : cs.incidentCells(e, 2))
                if(indexAttr[f].label == label)
                    if(normSum > 0) {
                        if(parm("Pin1 Sensitivity Average Method") == "Soft-max")
                          eD.pin1Sensitivity[label] = exp(eD.pin1SensitivityRaw[label]) / normSum;
                        else if(parm("Pin1 Sensitivity Average Method") == "Arithmetic Average")
                            eD.pin1Sensitivity[label] = eD.pin1SensitivityRaw[label] / normSum;
                    }

        }
    }

    // update miscellaneous cell chemicals attributes
    double avg_auxin_conc = 0; // needed for later debug print
    for(uint i = 0; i < cellAttr.size(); i++) {
        auto it = cellAttr.begin();
        advance(it, i);
        Tissue::CellData& cD = it->second;
        cD.auxinDecayRate = parm("Auxin Decay Rate").toDouble();
        if(cD.type == Tissue::QC) {
            cD.auxinProdRate = parm("Auxin QC Basal Production Rate").toDouble();
            for(CCIndex vn : csDual.neighbors(cD.dualVertex)) {
                int labeln = indexAttr[vn].label;
                Tissue::CellData& cDn = cellAttr[labeln];
                cDn.auxinProdRate = parm("Auxin SCN Basal Production Rate").toDouble();
            }
        }
        if(cD.type == Tissue::Substrate)
            cD.auxin = 0;
        // negative check of chemicals
        std::map<const char*, double*> chems = {
            {"Cell auxin", &(cD.auxin)}, {"Cell Aux1", &(cD.Aux1)}, {"Cell Pin1", &(cD.Pin1)}};
        for(auto chem : chems)
            if(*chem.second < 0)
                *chem.second = 0;
        avg_auxin_conc += cD.auxin/cD.area;
    }
    avg_auxin_conc /= cellAttr.size();

    // update miscellaneous edge chemicals (only on the walls)
    for(uint i = 0; i < tissueProcess->wallEdges.size(); i++) {
        auto p = tissueProcess->wallEdges.begin();
        advance(p, i);
        for(CCIndex e : p->second) {
            Tissue::EdgeData& eD = edgeAttr[e];
            // negative check of chemicals
            std::map<const char*, double*> chems = {
                {"Edge Intercellular Auxin", &(eD.intercellularAuxin)}};
            for(auto chem : chems)
                if(*chem.second < 0)
                    *chem.second = 0;
            for(auto& q : eD.Pin1)
                if(q.second < 0)
                    q.second = 0;
            for(auto& q : eD.Aux1)
                if(q.second < 0)
                    q.second = 0;
        }
    }

    // divide all debugs infos by the total vertices, so getting the mean
    for(auto p : debugs)
        p.second /= mesh->ccStructure("TissueDual").vertices().size();

    userTime += Dt;
    double normal = calcNorm();

    // Print some debug info
    if(++debug_step >=  parm("Debug Steps").toInt()) {
        if(parm("Verbose") == "True")
            mdxInfo << "Chemical Time: " << userTime <<" Norm: " << normal << " Average Auxin Conc. :" << avg_auxin_conc << endl;
        debug_step = 0;
    }  

    // Chemical convergence is ***only*** based on auxin flow
    if(normal <= convergeThresh) {
        if(parm("Verbose") == "True") {
            mdxInfo << "Chemical solver converged" <<  endl;
            return true;
        }
    } else
        return false;


    return false;
}

bool Chemicals::step() {
    return update();
}

void Tissue::Subdivide::splitCellUpdate(Dimension dim,
                                        const CCStructure& cs,
                                        const CCStructure::SplitStruct& ss,
                                        CCIndex otherP,
                                        CCIndex otherN,
                                        double interpPos, bool cellDivision) {

    if(dim == 1) {
        if(!edgeAttr or !vtxAttr)
            mdxInfo << "MechanicsSubdivide::splitCellUpdate Edge attribute not set" << endl;
        else {
            double prevLength = (*edgeAttr)[ss.parent].prevLength;
            (*edgeAttr)[ss.childP].prevLength = interpPos * prevLength;
            (*edgeAttr)[ss.childN].prevLength = (1.0 - interpPos) * prevLength;
            double restLength = (*edgeAttr)[ss.parent].restLength;
            (*edgeAttr)[ss.childP].restLength = interpPos * restLength;
            (*edgeAttr)[ss.childN].restLength = (1.0 - interpPos) * restLength;
            double intercellularAuxin = (*edgeAttr)[ss.parent].intercellularAuxin;
            (*edgeAttr)[ss.childP].intercellularAuxin = interpPos * intercellularAuxin;
            (*edgeAttr)[ss.childN].intercellularAuxin = (1.0 - interpPos) * intercellularAuxin;
            for(auto p : (*edgeAttr)[ss.parent].Pin1) {
                (*edgeAttr)[ss.childP].Pin1[p.first] = interpPos * p.second;
                (*edgeAttr)[ss.childN].Pin1[p.first] = (1.0 - interpPos) *  p.second;
                (*edgeAttr)[ss.childP].Aux1[p.first] = interpPos * p.second;
                (*edgeAttr)[ss.childN].Aux1[p.first] = (1.0 - interpPos) *  p.second;
                (*edgeAttr)[ss.childP].PINOID[p.first] = interpPos * p.second;
                (*edgeAttr)[ss.childN].PINOID[p.first] = (1.0 - interpPos) *  p.second;
                (*edgeAttr)[ss.childP].PP2A[p.first] = interpPos * p.second;
                (*edgeAttr)[ss.childN].PP2A[p.first] = (1.0 - interpPos) *  p.second;
            }
            if(cellDivision)
                (*vtxAttr)[ss.membrane].divisionPoint = true;
            auto Pbounds = cs.edgeBounds(ss.childP);
            auto Nbounds = cs.edgeBounds(ss.childN);

            // not sure what the crap is all of this dirichlet stuff
            Point3u dirichletFirst;
            Point3u dirichletSecond;
            Point3u dirichletNew = Point3u(0, 0, 0);

            if(Pbounds.first == ss.membrane)
                dirichletFirst = (*vtxAttr)[Pbounds.second].dirichlet;
            else
                dirichletFirst = (*vtxAttr)[Pbounds.first].dirichlet;
            if(Nbounds.first == ss.membrane)
                dirichletSecond = (*vtxAttr)[Nbounds.second].dirichlet;
            else
                dirichletSecond = (*vtxAttr)[Nbounds.first].dirichlet;

            if(dirichletFirst.x() && dirichletSecond.x())
                dirichletNew.x() = 1;
            if(dirichletFirst.y() && dirichletSecond.y())
                dirichletNew.y() = 1;
            if(dirichletFirst.z() && dirichletSecond.z())
                dirichletNew.z() = 1;
            (*vtxAttr)[ss.membrane].dirichlet = dirichletNew;
        }
    } else if(dim == 2) {
        if(!faceAttr)
            mdxInfo << "Tissue::Subdivide::splitCellUpdate Face attribute not set" << endl;
        else {
            // Update faces
            FaceData& pfD = (*faceAttr)[ss.childP];
            FaceData& nfD = (*faceAttr)[ss.childN];
            FaceData& rfD = (*faceAttr)[ss.parent];

            pfD = nfD = rfD; // just copy all fields (SAFE?)
            // assign new available labels to the daughters
            std::set<int> labels = getAllLabelsFaces(cs, *indexAttr);
            for(int i = 1; i < INT_MAX; i++)
                if(labels.find(i) == labels.end()) {
                    (*indexAttr)[ss.childP].label = i;
                    break;
                }
            labels.insert((*indexAttr)[ss.childP].label);
            for(int i = 1; i < INT_MAX; i++)
                if(labels.find(i) == labels.end()) {
                    (*indexAttr)[ss.childN].label = i;
                    break;
                }
            // What to do with PIN? For gradient it gets recalculated at every step.

            // set them as daugheters for later
            CellData& pCD = (*cellAttr)[(*indexAttr)[ss.parent].label];
            pCD.daughters =
                std::make_pair((*indexAttr)[ss.childP].label, (*indexAttr)[ss.childN].label);
        }
    }
}

// Initialize to grab subdivider
bool RootDivide::initialize(QWidget* parent) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("RootDivide::step No current mesh"));
    CCIndexDataAttr& indexAttr = mesh->indexAttr();

    // Call base initialize
    if(!CellDivision::initialize(parent))
        return false;

    // Setup subdivision objects
    subdiv = Splitter(true);
    subdiv.mdx = MDXSubdivideX(*mesh);
    subdiv.mechanics =
        Tissue::Subdivide(indexAttr,
                          mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData"),
                          mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData"),
                          mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData"),
                          mesh->attributes().attrMap<int, Tissue::CellData>("CellData"));

    return true;
}

// a copy of the original MDX class, because that one crashes
MDXSubdivideX::MDXSubdivideX(Mesh &_mesh)
{
  mesh = &_mesh;
  indexAttr = &mesh->indexAttr();

  /*for(const QString &signal : mesh->signalAttrList())
    attrVec.push_back(&mesh->signalAttr<double>(signal)); //this creates warnings
    */
}

void MDXSubdivideX::splitCellUpdate(Dimension dim, const CCStructure &cs,
             const CCStructure::SplitStruct& ss, CCIndex otherP, CCIndex otherN, double interpPos)
{
  // Return if not set up
  if(!indexAttr) {
    mdxInfo << "MDXSubdivide:splitCellUpdate Error, indexAttr not set" << endl;
    return;
  }
  // First do membrane part
  // Don't propagate MDX info for edges, it will generate tons of attributes that aren't used
  if(dim != 2 and !otherP.isPseudocell() and !otherN.isPseudocell()) {
    CCIndexData &dP = (*indexAttr)[otherP], &dN = (*indexAttr)[otherN], &dV = (*indexAttr)[ss.membrane];

    // Propagate the label
    if(dN.label > 0)
      dV.label = dN.label;
    else if(dP.label > 0)
      dV.label = dP.label;

    // Propagate selection
    if(dN.selected and dP.selected)
      dV.selected = true;
    else
      dV.selected = false;

    // Update signal maps for membrane
    /*for(auto attr : attrVec) {
        if(attr != 0)
            if((*attr)[otherP] != 0 && (*attr)[otherN] != 0)
                (*attr)[ss.membrane] = (1.0 - interpPos) * (*attr)[otherP] + interpPos * (*attr)[otherN];
    }*/
  }

  // Now do cell part

  // Don't propagate MDX info for edges, it will generate tons of attributes that aren't used
  if(dim != 1) {
    CCIndexData &eP = (*indexAttr)[ss.childP], &eN = (*indexAttr)[ss.childN], &ePar = (*indexAttr)[ss.parent];
    // Update signal maps for cell
    /*for(auto attr : attrVec)
      (*attr)[ss.childP] = (*attr)[ss.childN] = (*attr)[ss.parent];*/
    eP.label = eN.label = ePar.label;
    eP.selected = eN.selected = ePar.selected;
    eP.measure = (1. - interpPos) * ePar.measure;
    eN.measure = interpPos * ePar.measure;
  }

  if(nextSubdivide)
    nextSubdivide->splitCellUpdate(dim,cs,ss,otherP,otherN,interpPos);
}


void Splitter::splitCellUpdate(Dimension dim,
                                    const CCStructure& cs,
                                    const CCStructure::SplitStruct& ss,
                                    CCIndex otherP,
                                    CCIndex otherN,
                                    double interpPos) {

    mdx.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
    mechanics.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos, cellDivision);
}

// Run a step of cell division
bool RootDivide::step(double Dt) {
    if(parm("Cell Division Enabled") == "True")
        // Pass our subdivide
        return CellDivision::step(getMesh("Mesh 1"), &subdiv);
    else
        return false;
}

// copy here the modified MDXProcessCellDivide
// Divide this so it can be called on any tissue
bool CellDivision::initialize(QWidget* parent) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("CellTissueCellDivideOvile::initialize No current mesh"));

    Verbose = parm("Verbose") == "True";

    if(Verbose)
        mdxInfo << "CellDivision::initialize" << endl;

    processParms();

    if(!getProcess(parm("Root Process"), rootProcess))
        throw(QString("Root::initialize Cannot make root process"));
    if(!getProcess(parm("Remesh"), remeshProcess))
        throw(QString("Root::initialize Cannot make Remesh") + parm("Remesh"));
    if(!getProcess(parm("Triangulate Faces Process"), triangulateProcess))
        throw(QString("Root::initialize Cannot make Triangulate Faces Process") +
              parm("Triangulate Faces"));
    if(!getProcess(parm("ClearCells Process"), clearCellsProcess))
        throw(QString("Root::initialize Cannot make ClearCells") + parm("ClearCells Process"));
    if(!getProcess(parm("Tissue Process"), tissueProcess))
        throw(QString("Root::initialize Cannot make Tissue Process") + parm("Tissue Process"));


    // Just use the basic MDX data for subdivision
    subdiv = MDXSubdivideX(*mesh);

    return true;
}

bool CellDivision::processParms() {
    // Process parameters
    Verbose = stringToBool(parm("Verbose"));

    return true;
}

// Helper function: Ensure that v is not within distance mw of v1 or v2
double pushToWallMin(Point3d &v, Point3d v1, Point3d v2, double mw)
{
  Point3d v1v2 = v2 - v1;
  double v1v2l = norm(v1v2);
  v1v2 /= v1v2l;

  // len is the distance along the edge from v to v1
  double len = (v - v1) * v1v2, len1 = len;
  if(v1v2l <= mw * 2.0)
    // If the edge is too short, we put v at its midpoint
    len1 = (v1v2l/2.0);
  else if(len < mw)
    // If v is too close to v1, we put it at the minimum distance
    len1 = mw;
  else if(len > v1v2l - mw)
    // If v is too close to v2, we put it at the minimum distance
    len1 = v1v2l - mw;

  v = v1 + v1v2 * len1;
  return(fabs(len - len1));
}

// Structure representing a division wall in 2D
struct DivWall2d
{
  struct {
    Point3d pos;    // position of endpoint
    CCIndex vA, vB; // vertices of wall to be split
    double sfrac;   // fraction of position along wall from A to B
  } endpoints[2];
};

bool lineLineIntersectX(Point3d p1,  Point3d p2,  Point3d  q1,  Point3d q2,
                                                                    double &s, double &t, Point3d &intersect, double eps = 1e-5)
{
  intersect = Point3d(0,0,0);
  Point3d rp = p2 - p1;
  Point3d rq = q2 - q1;

  Point3d rxs = rp % rq;
  if(norm(rxs) >= eps) {
    Point3d qps = (q1-p1) % rq;
    Point3d qpr = (q1-p1) % rp;

    if(abs(rxs.x())>abs(rxs.y()) and abs(rxs.x())>abs(rxs.z())){
      s = qps.x()/rxs.x();
      t = qpr.x()/rxs.x();
    } else if(abs(rxs.y())>abs(rxs.z())){
      s = qps.y()/rxs.y();
      t = qpr.y()/rxs.y();
    } else {
      s = qps.z()/rxs.z();
      t = qpr.z()/rxs.z();
    }
    intersect = p1 + s*rp;
    return true;
  }

  return false;
}



/* Find cell division wall */
bool findCellDiv2dX(const CCStructure &cs, const CCIndexDataAttr &indexAttr, CCIndex cell,
                   const Cell2dDivideParms &divParms, Cell2dDivideParms::DivAlg divAlg,
                   const Point3d &divVector, DivWall2d &divWall)
{
    // Read parameters from Parms
    double CellWallMin = divParms.parm("Cell Wall Min").toDouble();
    double CellWallSample = divParms.parm("Cell Wall Sample").toDouble();
    double CellPinch = divParms.parm("Cell Pinch").toDouble();
    double CellMaxPinch = divParms.parm("Cell Max Pinch").toDouble();

    for(uint i = 0 ; i < 2 ; i++)
    divWall.endpoints[i].vA = divWall.endpoints[i].vB = CCIndex::UNDEF;

    // Check dimension
    Dimension dimension = cs.dimensionOf(cell);
    if(dimension != 2)
    throw(QString("findCellDiv2d: Can only split a two-dimensional cell!"));

    // iterate to get points on walls
    CCStructure::CellTuple ct(cs,cell);
    std::vector<Point3d> pos;
    CCIndex firstV = ct[0];
    do {
    pos.push_back(indexAttr[ct[0]].pos);
    ct.flip(0,1);
    } while(ct[0] != firstV);
    uint numVertices = pos.size();

    // get centroid and normal
    auto cdata = polygonCentroidData(pos);

    double shortestDist = HUGE_VAL;
    double smallestAngle = HUGE_VAL;

    // try walls from each position in turn
    for(uint i = 0 ; i < numVertices ; i++) {

      uint i1 = (i+1) % numVertices;
      Point3d u1 = pos[i], u2 = pos[i1], u1u2 = u2 - u1;
      double u1u2l = norm(u1u2);
      Point3d u1u2d = u1u2 / u1u2l;

      // Try a wall from points along the wall separated by CellWallSample
      double ulEnd = u1u2l - CellWallMin;
      for(double ul = CellWallMin ; ul <= ulEnd ; ul += CellWallSample) {
        Point3d u = u1 + ul * u1u2d, uc = cdata.centroid - u;
        Point3d nu = cdata.normal ^ uc;
        double centroidDist = norm(uc);

        // Now we check the upcoming walls for an intersection
        CCStructure::CellTuple ct1(ct);
        ct1.flip(0,1);
        for(uint j = i+1; j < numVertices; j++) {
          uint j1 = (j+1) % numVertices;
          Point3d v;
          double s;
          if(planeLineIntersect(u, nu, pos[j], pos[j1], s, v) and (s >= 0. && s <= 1.)) {
            if(divAlg == Cell2dDivideParms::SHORTEST_WALL_THROUGH_CENTROID) {
                double dist = norm(u - v);
                // Check that the line actually goes through the centroid
                // (needed if the cell is not convex)
                if(dist < centroidDist) continue;
                dist += pushToWallMin(v, pos[j], pos[j1], CellWallMin);
                if(dist < shortestDist) {
                  shortestDist = dist;
                  divWall.endpoints[0] = { u,  ct[0],  ct.other(0), norm(u - pos[i]) / norm(pos[i1] - pos[i]) };
                  divWall.endpoints[1] = { v, ct1[0], ct1.other(0), norm(v - pos[j]) / norm(pos[j1] - pos[j]) };
                }
            } else if (divAlg == Cell2dDivideParms::ASSIGNED_VECTOR_TRHOUGH_CENTROID) {
                double dist = norm(u - v);
                if(dist < centroidDist) continue;
                dist += pushToWallMin(v, pos[j], pos[j1], CellWallMin);
                Point3d vu = u   - v;
                Point3d vv = v   - u;
                double angleu = acos((vu * divVector) / (norm(vu) * norm(divVector)));
                double anglev = acos((vv * divVector) / (norm(vv) * norm(divVector)));

                if(angleu < smallestAngle ) {
                  smallestAngle = angleu;
                  shortestDist = dist;
                  divWall.endpoints[0] = { u,  ct[0],  ct.other(0), norm(u - pos[i]) / norm(pos[i1] - pos[i]) };
                  divWall.endpoints[1] = { v, ct1[0], ct1.other(0), norm(v - pos[j]) / norm(pos[j1] - pos[j]) };
                }  else if(anglev < smallestAngle ) {
                    smallestAngle = anglev;
                    shortestDist = dist;
                    divWall.endpoints[0] = { u,  ct[0],  ct.other(0), norm(u - pos[i]) / norm(pos[i1] - pos[i]) };
                    divWall.endpoints[1] = { v, ct1[0], ct1.other(0), norm(v - pos[j]) / norm(pos[j1] - pos[j]) };
                  }

            } else
                throw(QString("findCellDiv2dX: Unknown division algorithm"));
          }
          ct1.flip(0,1);
        }
      }
      ct.flip(0,1);
    }

    // If the CellWallMin is too high, maybe we can't find a wall
    if(divWall.endpoints[0].vA.isPseudocell() or divWall.endpoints[0].vB.isPseudocell() or
     divWall.endpoints[1].vA.isPseudocell() or divWall.endpoints[1].vB.isPseudocell()) {
    mdxInfo << "CellTissue::findCellDivX Unable to find divide wall for cell:"
            << to_string(cell).c_str() << ", label:" << indexAttr[cell].label
            << ", CellWallMin or CellWallSample too high?" <<  endl;
    return false;
    }

    // Now we pinch the walls.
    for(uint i = 0 ; i < 2 ; i++) {
    auto &ep = divWall.endpoints[i];
    CCIndex edge = cs.join(ep.vA,ep.vB);
    // Don't pinch a given cell wall if it is external
    if(!cs.onBorder(edge)) {
      Point3d centroidDir = cdata.centroid - ep.pos;
      double cdl = norm(centroidDir);
      centroidDir /= cdl;

      Point3d posA = indexAttr[ep.vA].pos, posB = indexAttr[ep.vB].pos;
      double distA = norm(ep.pos - posA), distB = norm(ep.pos - posB);
      double minD = std::min(distA , distB);
      double pinchAmt = std::min(minD , cdl) * CellPinch;
      pinchAmt = std::min(pinchAmt , CellMaxPinch);

      ep.pos += centroidDir * pinchAmt;
    }
    }

    return true;
}


CCIndex getClosestAvailableCutPoint(CCIndex v, Point3d& closest, const CCStructure &cs, const CCIndexDataAttr &indexAttr, CCIndex cell, double CellWallMin, double CellWallSample) {

    closest = Point3d(BIG_VAL, 0, 0);
    cout << "finding closest point to " << v << endl;
    Point3d v_pos = indexAttr[v].pos;
    CCStructure::CellTuple ct(cs,cell);
    CCIndex edge;
    CCIndex firstV = ct[0];
    do {
        Point3d u1 = indexAttr[ct[0]].pos;
        Point3d u2 = indexAttr[ct.other(0)].pos;
        Point3d u1u2 = u2 - u1;
        double u1u2l = norm(u1u2);
        Point3d u1u2d = u1u2 / u1u2l;
        double ulEnd = u1u2l - CellWallMin;
        for(double ul = CellWallMin ; ul <= ulEnd ; ul += CellWallSample) {
          Point3d u = u1 + ul * u1u2d;
          if(norm(u - v_pos) < norm(closest - v_pos)) {
              closest = u;
              edge = ct[1];
          }
        }
        ct.flip(0,1);
    } while(ct[0] != firstV);
    cout << "it's " << closest << " in " << edge << endl;
    return edge;
}

// Divide a two-dimensional cell using to the given division parameters and algorithm.
bool divideCell2dX(CCStructure &cs, CCIndexDataAttr &indexAttr, Tissue::VertexDataAttr &vMAttr, CCIndex cell,
                       const Cell2dDivideParms &divParms, Subdivide *sDiv,
                   Cell2dDivideParms::DivAlg divAlg = Cell2dDivideParms::SHORTEST_WALL_THROUGH_CENTROID,
                   const Point3d &divVector = Point3d(1., 0., 0.), std::set<CCIndex> divisionPoints = std::set<CCIndex>(), double maxJoiningDistance = 1)
{

  // Find out where the division will take place
  DivWall2d divWall;

  CCIndex ep[2];
  //if(!findCellDiv2dX(cs, indexAttr, cell, divParms, divAlg, divVector, divWall))
  //  return false;
  if(divAlg == Cell2dDivideParms::SHORTEST_WALL_THROUGH_CENTROID ) {
      if(!findCellDiv2dX(cs, indexAttr, cell, divParms, divAlg, divVector, divWall))
          return false;
  } else if (divAlg == Cell2dDivideParms::ASSIGNED_VECTOR_TRHOUGH_CENTROID) {
      if(!findCellDiv2dX(cs, indexAttr, cell, divParms, divAlg, divVector, divWall))
          return false;
  } else if     (divAlg == 2) {
      if(divisionPoints.empty())
          return(divideCell2dX(cs, indexAttr, vMAttr, cell, divParms, sDiv, Cell2dDivideParms::ASSIGNED_VECTOR_TRHOUGH_CENTROID, divVector));
      else {
          if(!findCellDiv2dX(cs, indexAttr, cell, divParms, Cell2dDivideParms::ASSIGNED_VECTOR_TRHOUGH_CENTROID, divVector, divWall))
              return false;
          if(divisionPoints.size() == 1) {
              CCIndex v = *(divisionPoints.begin());
              if(norm(divWall.endpoints[0].pos - indexAttr[v].pos) <  norm(divWall.endpoints[1].pos - indexAttr[v].pos)) {
                  if(norm(divWall.endpoints[0].pos - indexAttr[v].pos) < maxJoiningDistance)
                    ep[0] = v;
              } else {
                  if(norm(divWall.endpoints[1].pos - indexAttr[v].pos) < maxJoiningDistance)
                    ep[1] = v;
              }
              vMAttr[v].divisionPoint = false;
          } else {
              ep[0] = *(divisionPoints.begin());
              for(CCIndex v : divisionPoints)
                  if(norm(divWall.endpoints[0].pos - indexAttr[v].pos) <  norm(divWall.endpoints[0].pos - indexAttr[ep[0]].pos)) {
                       ep[0] = v;
                       vMAttr[v].divisionPoint = false;
                  }
              if(norm(divWall.endpoints[0].pos - indexAttr[ep[0]].pos) > maxJoiningDistance)
                  ep[0] = CCIndex::UNDEF;
              ep[1] = *(divisionPoints.begin());
              for(CCIndex v : divisionPoints)
                  if(norm(divWall.endpoints[1].pos - indexAttr[v].pos) <  norm(divWall.endpoints[1].pos - indexAttr[ep[1]].pos)) {
                       ep[1] = v;
                       vMAttr[v].divisionPoint = false;
                  }
              if(norm(divWall.endpoints[1].pos - indexAttr[ep[1]].pos) > maxJoiningDistance)
                  ep[1] = CCIndex::UNDEF;
          }

      }
  } else
      throw(QString("divideCell2dX: Unknown division algorithm"));

  // Divide the cell walls
  if(divAlg == Cell2dDivideParms::SHORTEST_WALL_THROUGH_CENTROID || divAlg == Cell2dDivideParms::ASSIGNED_VECTOR_TRHOUGH_CENTROID) {
      for(int i = 0 ; i < 2 ; i++) {
        CCIndex edge = edgeBetween(cs, divWall.endpoints[i].vA, divWall.endpoints[i].vB);

        CCStructure::SplitStruct ss(edge);
        CCIndexFactory.fillSplitStruct(ss);
        ep[i] = ss.membrane;

        bool neg = cs.ro(cell, edge) == ccf::NEG;
        cs.splitCell(ss);
        // Update the position ourselves since the subdivider may need it
        indexAttr[ep[i]].pos = divWall.endpoints[i].pos;
        if(sDiv) {
          sDiv->splitCellUpdate(1, cs, ss, divWall.endpoints[i].vA, divWall.endpoints[i].vB,
                                neg ? divWall.endpoints[i].sfrac : (1.0 - divWall.endpoints[i].sfrac));
        }
      }
  } else if(divAlg == 2) {
      for(int i = 0 ; i < 2 ; i++) {
        if(!ep[i].isPseudocell()) continue;

        CCIndex edge = edgeBetween(cs, divWall.endpoints[i].vA, divWall.endpoints[i].vB);

        CCStructure::SplitStruct ss(edge);
        CCIndexFactory.fillSplitStruct(ss);
        ep[i] = ss.membrane;

        bool neg = cs.ro(cell, edge) == ccf::NEG;
        cs.splitCell(ss);
        // Update the position ourselves since the subdivider may need it
        indexAttr[ep[i]].pos = divWall.endpoints[i].pos;
        if(sDiv) {
          sDiv->splitCellUpdate(1, cs, ss, divWall.endpoints[i].vA, divWall.endpoints[i].vB,
                                neg ? divWall.endpoints[i].sfrac : (1.0 - divWall.endpoints[i].sfrac));
        }
      }
  }

  // Divide the cell itself
  CCStructure::SplitStruct ss(cell);
  CCIndexFactory.fillSplitStruct(ss);
  cs.splitCell(ss, +ep[0] -ep[1]);

  if(sDiv) {
    updateFaceGeometry(cs, indexAttr, ss.childP);
    updateFaceGeometry(cs, indexAttr, ss.childN);
    double s = indexAttr[ss.childP].measure / (indexAttr[ss.childP].measure + indexAttr[ss.childN].measure);
    sDiv->splitCellUpdate(2, cs, ss, CCIndex(), CCIndex(), s);
  }

  return true;
}

// Run a step of cell division
bool CellDivision::step(Mesh* mesh, Subdivide* subdiv) {
    if(!mesh)
        throw(QString("CellDivision.step: Mesh not initialized"));

    CCStructure& cs = mesh->ccStructure(mesh->ccName());
    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    Tissue::FaceDataAttr& faceAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");

    bool manualCellDivision = parm("Manual Cell Division Enabled") == "True";
    double divisionMeristemSize = parm("Division Meristem Size").toDouble();
    double divisionMaxTime = parm("Max Division Time").toDouble();
    double divisionMinTime = parm("Min Division Time").toDouble();
    double divisionProbHalfSize = parm("Division half-probability by Cell Size Ratio").toDouble();
    double divisionProbHalfInhibitor = parm("Division half-probability by Inhibitor").toDouble();
    double divisionPromoterLevel = parm("Division Promoter Level").toDouble();
    bool divisionControl = parm("Division Control") == "True";
    bool ignoreCellType = parm("Ignore Cell Type") == "True";

    // find the QC so we can print the distance (for plotting)
    Point3d  QCcm;
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        if(cD.type == Tissue::QC)
            QCcm += cD.centroid;
    }
    QCcm /= 2;
    // check if there are any cells to divide (depending on are or other clues)
    std::vector<Tissue::CellData> cDs;
    bool trigger_division = false;
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        cD.divProb = 0;
        cD.divisionAllowed = false;
        cD.divProb = (2 / (1 + exp(-divisionProbHalfInhibitor * cD.divInhibitor*100/cD.area)) - 1) * divisionMaxTime + divisionMinTime;
        if(divisionControl && rootProcess->userTime > 24) {
            if(cD.area/cD.cellMaxArea > divisionProbHalfSize)
                if(norm(cD.centroid - QCcm) < divisionMeristemSize)
                    if(cD.divPromoter/cD.area > divisionPromoterLevel) {
                        if(cD.lastDivision > cD.divProb )
                            cD.divisionAllowed = true;
                    }
        } else
            cD.divisionAllowed = true;

        if((manualCellDivision && cD.selected) ||
            (cD.divisionAllowed == true &&
             cD.area > cD.cellMaxArea &&
             cD.lastDivision > 1)) {
            if(Verbose) {
                mdxInfo << "CellDivision.step: Cell division triggered by " << cD.label << " of size " << cD.area
                        << " of type " << Tissue::ToString(cD.type) << " at position " << cD.centroid << " distance from QC " << cD.centroid.y() - QCcm.y()
                        <<   " bigger than " << cD.cellMaxArea << " last division time: " << cD.lastDivision
                        << " division inhibitor: " << cD.divInhibitor/cD.area  << " division promoter: " << cD.divPromoter/cD.area  << " division prob: " << cD.divProb  << endl;


                //cout << random << " " << divProbAuxin << " " << divProbSize << " " << divProbInhibitor << " " << (divProbSize + divProbAuxin *divProbInhibitor* divProbSize)*0.5 << endl;
                //cout << (random < RAND_MAX * (divProbSize + divProbAuxin * divProbSize * divProbInhibitor)*0.5 * Dt) << endl;

            }
            // Skip division if division algoritm MF depending but cell has no polarity
            if(cD.divAlg != 0 && norm(cD.a1) < parm("Minimum Polarity Vector Norm").toDouble()) {
                if(Verbose)
                    mdxInfo << "CellDivision.step: WARNING: Non-polar division prevented: " << cD.label << " division vector norm: " << norm(cD.a1) << endl;
                 continue;
            }
            trigger_division = true;
            cDs.push_back(cD);
            break;  ///// Known bug: dividing two adjacient cells at the same time provokes the Undefined type cell bug
        }
    }

    // no cell to divide, return
    if(!trigger_division)
        return false;

    // there is at least one cell divide, turn the cells into single faces and
    // find again which ones to divide
    clearCellsProcess->clearCell(cDs[0].label);
    updateGeometry(cs, indexAttr);
    std::vector<CCIndex> D;
    for(CCIndex f : cs.faces()) {
        Tissue::CellData& cD = cellAttr[indexAttr[f].label];
        if((find(cDs.begin(), cDs.end(), cD) != cDs.end() ))
            D.push_back(f);
    }

    // something goes, wrong, we cannot find the cells to divide
    if(D.empty())
        throw(QString("CellDivision.step: empty cells division set"));

    // divide the cells
    std::map<int, std::pair<int, int>> daughters;        
    forall(const CCIndex& f, D) {
        int label = indexAttr[f].label;
        Tissue::CellData& cD = cellAttr[label];
        uint DivAlg = parm("Division Algorithm").toInt();
        if(cD.divAlg > -1)
            DivAlg = cD.divAlg;
        if(Verbose)
            mdxInfo << "CellDivision.step: Dividing cell (" << label << "), of type " << Tissue::ToString(cD.type) << " with area " << indexAttr[f].measure << " using algorithm: " << DivAlg <<  endl;        // Simple division, shortest edge through centroid
        // Division along closest walls
        if(DivAlg == 0) {
               if(!divideCell2dX(cs, indexAttr, vMAttr, f, *this, subdiv))
                    throw(QString("CellDivision.step: Failed division for: " +  label));
        // Division along vector through centroid
        } else if(DivAlg == 1){
            if(Verbose)
                mdxInfo << "CellDivision.step: Dividing cell along : " << cD.divVector << endl;
            if(!divideCell2dX(cs, indexAttr, vMAttr, f, *this, subdiv, ASSIGNED_VECTOR_TRHOUGH_CENTROID, cD.divVector)){
                mdxInfo << "CellDivision.step: Failed division for: " +  label << endl;
                return false;
            }
        // Division point assisted division
        } else if(DivAlg == 2) {
            if(Verbose)
                mdxInfo << "CellDivision.step: Dividing cell assisted by division points : " << endl;
            std::set<CCIndex> divisionPoints;
            for(CCIndex v : cD.perimeterVertices)
                if(vMAttr[v].divisionPoint)
                    divisionPoints.insert(v);
            if(!divideCell2dX(cs, indexAttr, vMAttr, f, *this, subdiv, (mdx::Cell2dDivideParms::DivAlg)2, cD.divVector,
                              divisionPoints, parm("Max Joining Distance").toDouble())){
                throw(QString("CellDivision.step: Failed division for: " +  label));
            }

        }
          else
            throw(QString("CellDivision.step: Unknown Division Algorithm") + DivAlg);

        if(&cellAttr[label] != nullptr)
            daughters[label] = cellAttr[label].daughters;
        else
            throw(QString("CellDivision.step: null cell" +  label));
    }

    // re-triangulate the cell, recreate the dual tissue and soft remesh
    indexAttr[*D.begin()].selected = true;
    triangulateProcess->step();
    tissueProcess->initialize();
    remeshProcess->step(true, false, true);

    // update the daughter cells
    std::map<Tissue::CellType, int> maxAreas;
    for (int fooInt = Tissue::Undefined; fooInt != Tissue::Source; fooInt++ )
        maxAreas[(Tissue::CellType)fooInt] =
                rootProcess->setGlobalAttrProcess->parm(QString(Tissue::ToString((Tissue::CellType)fooInt)) + " Max area").toInt();
    for(Tissue::CellData cD : cDs)
        if(daughters[cD.label].first > 0 && daughters[cD.label].second > 0)
            cD.division(cs, cellAttr, faceAttr, edgeAttr,
                        cellAttr[daughters[cD.label].first], cellAttr[daughters[cD.label].second], maxAreas, ignoreCellType);

    // Update mesh points, edges, surfaces
    mesh->updateAll();

    return true;
}


void write_debug_header(std::ofstream &output_file) {
    output_file << "Steps"
                << ","
                << "User Time"
                << ","
                << "Real Time"
                << ","
                << "Growth Rate"
                << ","
                << "QC position"
                << ","
                << "Average Pressure"
                << ","
                << "Average Sigma"
                << ","
                << "Average Cell Growth"
                << ","
                << "Average Vertex Velocity"
                << ","
                << "Average Auxin Cytoplasm"
                << ","
                << "Std Auxin Cytoplasm"
                << ","
                << "Average Auxin Intercellular"
                << ","
                << "Average PIN Cytoplasm"
                << ","
                << "Average PIN Membrane"
                << ","
                << "Average AUX1 Cytoplasm"
                << ","
                << "Average Auxin Diffusion"
                << ","
                << "Average Auxin Export"
                << ","
                << "Average PIN1 Expression"
                << ","
                << "Average PIN1 Trafficking"
                << ","
                << "Anisotropy degree"
                << endl
                << flush;
}

// Initialize the main solver process
bool Root::initialize(QWidget* parent) {
    srand (1983);

    mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("Root::initialize No current mesh"));
    //CCStructure& cs = mesh->ccStructure("Tissue");
    //CCIndexDataAttr& indexAttr = mesh->indexAttr();

    this->widget_parent = parent;

    if(!getProcess(parm("Tissue Process"), tissueProcess))
        throw(QString("Root::initialize Cannot make Tissue Process:") + parm("Tissue Process"));
    if(!getProcess(parm("Mechanics Process"), mechanicsProcess))
        throw(QString("Root::initialize Cannot make Mechanics:") + parm("Mechanics"));
    if(!getProcess(parm("Mechanical Growth Process"), mechanicalGrowthProcess))
        throw(QString("Root::initialize Cannot make Mechanical Growth Process:") +
              parm("Mechanical Growth Process"));
    if(!getProcess(parm("Chemicals Process"), chemicalsProcess))
        throw(QString("Root::initialize Cannot make Chemicals Process:") +
              parm("Chemicals Process"));
    if(!getProcess(parm("Divide Process"), divideProcess))
        throw(QString("Root::initialize Cannot make Divide Process:") + parm("Divide Process"));
    if(!getProcess(parm("Delete Cell Process"), deleteProcess))
        throw(QString("Root::initialize Cannot make Delete Cell Process:") + parm("Delete Cell Process"));
    if(!getProcess(parm("Execute Test Process"), executeTestProcess))
        throw(QString("Root::initialize Cannot make Execute Test Process:") + parm("Execute Test Process"));
    if(!getProcess(parm("Set Global Attr Process"), setGlobalAttrProcess))
        throw(QString("Root::initialize Cannot make Set Global Attr Process:") +
              parm("Set Global Attr Process"));
    if(!getProcess(parm("Triangulate Faces"), triangulateProcess))
        throw(QString("Root::initialize Cannot make Triangulate Faces") +
              parm("Triangulate Faces"));
    if(!getProcess(parm("SplitEdges"), splitEdgesProcess))
        throw(QString("Root::initialize Cannot make SplitEdges") + parm("SplitEdges"));
    if(!getProcess(parm("SaveView"), saveViewProcess))
        throw(QString("Root::initialize Cannot make SaveView") + parm("SaveView"));
    if(!getProcess(parm("SaveMesh"), saveMeshProcess))
        throw(QString("Root::initialize Cannot make SaveMesh") + parm("SaveMesh"));
    if(!getProcess(parm("Remesh"), remeshProcess))
        throw(QString("Root::initialize Cannot make Remesh") + parm("Remesh"));

    debugging = parm("Debug") == "True";
    maxMechanicsIter =  parm("Max Mechanical Iterations").toInt();
    maxChemicalIter =  parm("Max Chemical Iterations").toInt();
    mechanicsEnabled =  parm("Mechanics Enabled") == "True";
    chemicalsEnabled = parm("Chemicals Enabled") == "True";
    growthEnabled =  parm("Growth Enabled") == "True";
    divisionEnabled = parm("Cell Division Enabled") == "True";

    // processes set
    processes.clear();
    processes.push_back(divideProcess);
    processes.push_back(mechanicsProcess);
    processes.push_back(mechanicalGrowthProcess);
    processes.push_back(chemicalsProcess);
    processes.push_back(setGlobalAttrProcess);
    // initialize the remesh process
    remeshProcess->initialize(parent);
    remeshProcess->setProcesses(processes);

    if(!process_started) {
        mdxInfo << "Initializing Root..." << endl;
        // remesh if necessary
        if(parm("Remesh at start") != "None")
            remeshProcess->step(true, false, parm("Remesh at start") == "Soft");
        // give it a grow
        mechanicalGrowthProcess->initialize(parent);
        mechanicalGrowthProcess->step(mechanicsProcess->Dt);
        // to avoid problems with the strain constraint
        mechanicsProcess->initialize(parent);
        mechanicsProcess->PBDProcess->rewind(parent);
        // init the debug file
        if(parm("Debug") == "True") {
            if(!output_file.is_open()) {
                output_file.open(parm("Debug File").toStdString());
                write_debug_header(output_file);
            }
        }
        // create the snapshots dir
        if(parm("Snapshots Timer").toInt() > 0) {
            if(system("mkdir \"snapshots\"") == -1)
                mdxInfo << "Error creating snapshot directory" << endl;
            /*if(system("mkdir \"views\"") == -1)
                mdxInfo << "Error creating view directory" << endl;*/
            if(parm("Snapshots Directory") == "")
                snapshotDir = "snapshots/" + currentDateTime() + "/";
            else
                snapshotDir = "snapshots/" + parm("Snapshots Directory").toStdString() + "/";
            //string viewFile = "views/" + currentDateTime() + ".mdxv";
            //saveView(QString::fromStdString(viewFile));
            if(system((string("mkdir ") + snapshotDir).c_str()) == -1)
                throw(QString("Error creating snapshot directory"));
        }
        begin_clock = clock();
        mdxInfo << "Root initialization completed." << endl;
        process_started = true;
    }

    mesh->updateAll();
    return true;
}


bool Root::rewind(QWidget* parent) {
    // To rewind, we'll reload the mesh
    mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("Root::rewind No current mesh, cannot rewind"));
    process_started = false;
    if(parm("Debug") == "True") {
        if(!output_file.is_open()) {
            output_file.open(parm("Debug File").toStdString());
            write_debug_header(output_file);
            output_file.close();
        } else {
            output_file.close();
            output_file.open(parm("Debug File").toStdString());
            write_debug_header(output_file);
        }

    }
    userTime = 0;
    stepCount = 0;
    screenShotCount = 0;
    mechanicsProcess->rewind(parent);
    chemicalsProcess->rewind(parent);
    executeTestProcess->rewind(parent);
    mechanicsProcessConverged = false;
    MeshLoad meshLoad(*this);
    meshLoad.setParm("File Name", mesh->file());
    return meshLoad.run();
}

// When calling --run on the console this does not work
// bool Root::run() { mdxInfo << "Running Root" << endl; return step(); }
bool Root::step() {

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
        throw(QString("Root::run Error, no cell complex selected"));
    if(ccName != "Tissue")
        throw(QString("Root::run Error, please tun the model on the Tissue Complex"));
    CCStructure& cs = mesh->ccStructure("Tissue");
    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");

    stepCount++;

    std::set<int> bad_cells;
    for(auto c : cellAttr)
        if(c.second.label == -1)
            bad_cells.insert(c.first);
    for(int index : bad_cells) {
        mdxInfo << "Bad cell found " << endl;
        cellAttr.erase(index);
    }

    // Update the mechanicals
    int i = 0;
    do{
        if(mechanicsEnabled)
            mechanicsProcessConverged= mechanicsProcess->step();
        userTime += mechanicsProcess->Dt;
        // Update the chemicals
        if(chemicalsEnabled)
            for(int j = 0; j < maxChemicalIter; j++)
                if(chemicalsProcess->update())
                    break;
        // Grow the tissue
        if(growthEnabled)
            mechanicalGrowthProcess->step(mechanicsProcess->Dt);
        // Perform cell division
        if(divisionEnabled)
            //while(divideProcess->step(mechanicsProcess->Dt)); // multiple divisions on the same step, usually crashes
            divideProcess->step(mechanicsProcess->Dt);
        // Reset attributes ?
        //setGlobalAttrProcess->step();
        // Update the tissue
        tissueProcess->step(mechanicsProcess->Dt);

        i++;
    } while(!mechanicsProcessConverged && i < maxMechanicsIter);

    // Execute tests
    executeTestProcess->step(stepCount);

    // Find the center of QC
    Point3d  QCcm;
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        if(cD.type == Tissue::QC)
            QCcm += cD.centroid;
    }
    QCcm /= 2;

    // Find the center of vascular initials or substrate
    Point3d  VIcm, SUBcm;
    int VIcount = 0;
    int SUBcount = 0;
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        if(cD.type == Tissue::VascularInitial) {
            VIcm += cD.centroid;
            VIcount++;
        }
        if(cD.type == Tissue::Substrate) {
            SUBcm += cD.centroid;
            SUBcount++;
        }
    }
    VIcm /= VIcount;
    SUBcm /= SUBcount;

    // Fix the camera on the Vascular initials
    QStringList list_QC = parm("Frame fixed on QC").split(QRegExp(","));
    if(list_QC.size() != 2)
         throw(QString("Frame fixed on QC should be two numbers separated by comma"));
    double x_QC = list_QC[0].toInt();
    double y_QC = list_QC[1].toInt();
    if(x_QC != 0 || y_QC != 0) {
        currentStack()->getFrame().setTranslation(-VIcm.x()-x_QC, -VIcm.y()-y_QC, 0);
        updateState();
        updateViewer();
    }

    // Fix the camera on the Substrate
    QStringList list_substrate = parm("Frame fixed on Substrate").split(QRegExp(","));
    if(list_substrate.size() != 2)
         throw(QString("Frame fixed on Substrate should be two numbers separated by comma"));
    double x_substrate = list_substrate[0].toInt();
    double y_substrate = list_substrate[1].toInt();
    if(x_substrate != 0 || y_substrate != 0) {
        currentStack()->getFrame().setTranslation(-SUBcm.x()-x_substrate, -SUBcm.y()-y_substrate, 0);
        updateState();
        updateViewer();
    }

    // Check whether we should remesh
    if(parm("Remesh during execution") == "True")
        remeshProcess->check();

    // Update the mesh
    if(stepCount % parm("Mesh Update Timer").toInt() == 0)
        mesh->updateAll();

    // Save a snapshot of the simulation
    if(parm("Snapshots Timer").toInt() > 0 && stepCount % parm("Snapshots Timer").toInt()  == 0){
        mdxInfo << "Let's take a snapshot" << endl;
        std::set<QString> signals_set = {
                                         "Chems: Division Inhibitor by Area",
                                         "Chems: Division Promoter by Area",
                                         "Chems: Division Probability",
                                         "Division Count",
                                         "Mechs: Growth Rate",
                                         "Chems: Auxin By Area"
                                        };
        for(QString signalName: signals_set) {
            mesh->updateProperties("Tissue");
            mesh->drawParms("Tissue").setGroupVisible("Faces", true);
            mesh->drawParms("Tissue").setRenderChoice("Faces", signalName);
            // move unwanted visual back
            if(signalName == QString("Chems: Auxin By Area") ||
                    signalName == QString("Chems: Division Promoter by Area") ||
                    signalName == QString("Chems: Division Inhibitor by Area") ||
                    signalName == QString("Chems: Division Time") ||
                    signalName == QString("Chems: Division Probability") ||
                    signalName == QString("Division Count")) {
                for(uint i = 0; i < cellAttr.size(); i++) {
                    auto it = cellAttr.begin();
                    advance(it, i);
                    Tissue::CellData& cD = it->second;
                    indexAttr[cD.PDGmax_e].pos[2] -= 10;
                    indexAttr[cD.PDGmax_v1].pos[2] -= 10;
                    indexAttr[cD.PDGmax_v2].pos[2] -= 10;
                    indexAttr[cD.PDGmin_e].pos[2] -= 10;
                    indexAttr[cD.PDGmin_v1].pos[2] -= 10;
                    indexAttr[cD.PDGmin_v2].pos[2] -= 10;
                }
            }
            if(signalName == QString("Mechs: Growth Rate") ||
                    signalName == QString("Chems: Division Promoter by Area") ||
                    signalName == QString("Chems: Division Inhibitor by Area") ||
                    signalName == QString("Chems: Division Time") ||
                    signalName == QString("Chems: Division Probability") ||
                    signalName == QString("Division Count")) {
                for(uint i = 0; i < cellAttr.size(); i++) {
                    auto it = cellAttr.begin();
                    advance(it, i);
                    Tissue::CellData& cD = it->second;
                    indexAttr[cD.auxinFlux_v1].pos[2] -= 10;
                    indexAttr[cD.auxinFlux_v2].pos[2] -= 10;
                    indexAttr[cD.auxinFlux_v3].pos[2] -= 10;
                    indexAttr[cD.auxinFlux_v4].pos[2] -= 10;
                    indexAttr[cD.auxinFlux_e].pos[2] -= 10;
                    indexAttr[cD.auxinFlux_el].pos[2] -= 10;
                    indexAttr[cD.auxinFlux_er].pos[2] -= 10;
                    indexAttr[cD.auxinFlux_eb].pos[2] -= 10;
                    indexAttr[cD.auxinFlux_f].pos[2] -= 10;
                }
                for(CCIndex e : cs.edges()) {
                    Tissue::EdgeData& eD = edgeAttr[e];
                    indexAttr[eD.pin_v1_1].pos[2] -= 10;
                    indexAttr[eD.pin_v2_1].pos[2] -= 10;
                    indexAttr[eD.pin_v3_1].pos[2] -= 10;
                    indexAttr[eD.pin_v4_1].pos[2] -= 10;
                    indexAttr[eD.pin_e1_1].pos[2] -= 10;
                    indexAttr[eD.pin_e2_1].pos[2] -= 10;
                    indexAttr[eD.pin_e3_1].pos[2] -= 10;
                    indexAttr[eD.pin_e4_1].pos[2] -= 10;
                    indexAttr[eD.pin_f_1].pos[2] -= 10;
                    indexAttr[eD.pin_v1_2].pos[2] -= 10;
                    indexAttr[eD.pin_v2_2].pos[2] -= 10;
                    indexAttr[eD.pin_v3_2].pos[2] -= 10;
                    indexAttr[eD.pin_v4_2].pos[2] -= 10;
                    indexAttr[eD.pin_e1_2].pos[2] -= 10;
                    indexAttr[eD.pin_e2_2].pos[2] -= 10;
                    indexAttr[eD.pin_e3_2].pos[2] -= 10;
                    indexAttr[eD.pin_e4_2].pos[2] -= 10;
                    indexAttr[eD.pin_f_2].pos[2] -= 10;
                }
            }

            mesh->setSignal(signalName);
            if(signalName == QString("Chems: Division Inhibitor by Area")) {
                double divInhibitor = divideProcess->parm("Division half-probability by Inhibitor").toDouble();
                if(divInhibitor > 1)divInhibitor=1;
                if(divInhibitor < 0)divInhibitor=0;
                if(divInhibitor > 1)
                    divInhibitor=10;
                else
                    divInhibitor=1;
                mesh->setSignalBounds(Point2d(0, divInhibitor));
            }
            mesh->updateAll();
            QString fileName = QString::fromStdString(snapshotDir) + QString("Root-%1-%2.png").arg(signalName).arg(screenShotCount, 4, 10, QChar('0'));
            takeSnapshot(fileName, 1, 645*4, 780*4, 10, true); // cluster?
            //takeSnapshot(fileName, 1, 2490*2, 1310*2, 100, true); // lab PC
            // restore unwanted visual forward
            if(signalName == QString("Chems: Auxin By Area") ||
                    signalName == QString("Chems: Division Promoter by Area") ||
                    signalName == QString("Chems: Division Inhibitor by Area") ||
                    signalName == QString("Chems: Division Time") ||
                    signalName == QString("Chems: Division Probability") ||
                    signalName == QString("Division Count")) {
                for(uint i = 0; i < cellAttr.size(); i++) {
                    auto it = cellAttr.begin();
                    advance(it, i);
                    Tissue::CellData& cD = it->second;
                    indexAttr[cD.PDGmax_e].pos[2] += 10;
                    indexAttr[cD.PDGmax_v1].pos[2] += 10;
                    indexAttr[cD.PDGmax_v2].pos[2] += 10;
                    indexAttr[cD.PDGmin_e].pos[2] += 10;
                    indexAttr[cD.PDGmin_v1].pos[2] += 10;
                    indexAttr[cD.PDGmin_v2].pos[2] += 10;
                }
            }
            if(signalName == QString("Mechs: Growth Rate") ||
                    signalName == QString("Chems: Division Promoter by Area") ||
                    signalName == QString("Chems: Division Inhibitor by Area") ||
                    signalName == QString("Chems: Division Time") ||
                    signalName == QString("Chems: Division Probability") ||
                    signalName == QString("Division Count")) {
                for(uint i = 0; i < cellAttr.size(); i++) {
                    auto it = cellAttr.begin();
                    advance(it, i);
                    Tissue::CellData& cD = it->second;
                    indexAttr[cD.auxinFlux_v1].pos[2] += 10;
                    indexAttr[cD.auxinFlux_v2].pos[2] += 10;
                    indexAttr[cD.auxinFlux_v3].pos[2] += 10;
                    indexAttr[cD.auxinFlux_v4].pos[2] += 10;
                    indexAttr[cD.auxinFlux_e].pos[2] += 10;
                    indexAttr[cD.auxinFlux_el].pos[2] += 10;
                    indexAttr[cD.auxinFlux_er].pos[2] += 10;
                    indexAttr[cD.auxinFlux_eb].pos[2] += 10;
                    indexAttr[cD.auxinFlux_f].pos[2] += 10;
                }
                for(CCIndex e : cs.edges()) {
                    Tissue::EdgeData& eD = edgeAttr[e];
                    indexAttr[eD.pin_v1_1].pos[2] += 10;
                    indexAttr[eD.pin_v2_1].pos[2] += 10;
                    indexAttr[eD.pin_v3_1].pos[2] += 10;
                    indexAttr[eD.pin_v4_1].pos[2] += 10;
                    indexAttr[eD.pin_e1_1].pos[2] += 10;
                    indexAttr[eD.pin_e2_1].pos[2] += 10;
                    indexAttr[eD.pin_e3_1].pos[2] += 10;
                    indexAttr[eD.pin_e4_1].pos[2] += 10;
                    indexAttr[eD.pin_f_1].pos[2] += 10;
                    indexAttr[eD.pin_v1_2].pos[2] += 10;
                    indexAttr[eD.pin_v2_2].pos[2] += 10;
                    indexAttr[eD.pin_v3_2].pos[2] += 10;
                    indexAttr[eD.pin_v4_2].pos[2] += 10;
                    indexAttr[eD.pin_e1_2].pos[2] += 10;
                    indexAttr[eD.pin_e2_2].pos[2] += 10;
                    indexAttr[eD.pin_e3_2].pos[2] += 10;
                    indexAttr[eD.pin_e4_2].pos[2] += 10;
                    indexAttr[eD.pin_f_2].pos[2] += 10;
                }
            }

        }
        screenShotCount++;
    }

    // Debugs?
    if(debugging) {
        double rootGR = 0;
        for(double gr : mechanicsProcess->growthRatesVector)
            rootGR += gr;
        rootGR /= mechanicsProcess->growthRatesVector.size();

        double avg_pressure = 0;
        double avg_sigma = 0;
        for(CCIndex e : cs.edges()) {
            Tissue::EdgeData eD = edgeAttr[e];
            avg_pressure += norm(eD.pressureForce);
            avg_sigma += norm(eD.sigmaForce);

        }
        avg_pressure *= 1. / cs.edges().size();
        avg_sigma *= 1. / cs.edges().size();


        double avg_velocity = 0;
        for(auto v : cs.vertices()) {
            Tissue::VertexData& vD = vMAttr[v];
            avg_velocity += norm(vD.velocity);
        }
        avg_velocity *= 1. / cs.vertices().size();

        double total_auxin = 0;
        double avg_auxin = 0, avg_auxinInter = 0, avg_Pin1Cyt = 0, avg_Pin1Mem = 0, avg_Aux1Mem = 0, avg_cell_growth = 0;
        for(auto c : cellAttr) {
            Tissue::CellData& cD = cellAttr[c.first];
            total_auxin += cD.auxin;
            avg_auxin += cD.auxin / cD.area;
            avg_Pin1Cyt += cD.Pin1 / cD.area;
            avg_cell_growth += cD.growthRate;
            for(CCIndex e : cD.perimeterEdges) {
                Tissue::EdgeData eD = edgeAttr[e];
                avg_Pin1Mem += eD.Pin1[cD.label] / eD.length;
                avg_Aux1Mem += cD.Aux1 / eD.length;
                avg_auxinInter += eD.intercellularAuxin / 2;
            }
        }
        avg_auxin *= 1. / cellAttr.size();
        avg_Pin1Cyt *= 1. / cellAttr.size();
        avg_Pin1Mem *= 1. / cs.edges().size();
        avg_Aux1Mem *= 1. / cs.edges().size();
        avg_auxinInter *= 1. / cs.edges().size();
        avg_cell_growth *= 1. / cellAttr.size();

        double std_auxin = 0;
        for(auto c : cellAttr) {
            Tissue::CellData& cD = cellAttr[c.first];
            std_auxin += ((cD.auxin / cD.area) - avg_auxin) * ((cD.auxin / cD.area) - avg_auxin);
        }
        std_auxin = sqrt(std_auxin / cellAttr.size());


        for(CCIndex e : cs.edges()) {
            Tissue::EdgeData eD = edgeAttr[e];
            if(eD.type == Tissue::Wall)
                total_auxin += eD.intercellularAuxin;
        }

        double anisotropy_degree = 0;
        int anisotropy_count = 0;
        for(auto c : cellAttr) {
            Tissue::CellData& cD = cellAttr[c.first];
            if(cD.mfRORate > 0 && cD.centroid[1] < 0) {
                anisotropy_degree += norm(cD.axisMax) / norm(cD.axisMin);
                anisotropy_count++;
            }
        }
        anisotropy_degree /= anisotropy_count;

        output_file << stepCount << "," << userTime << "," << mechanicsProcess->realTime  << ","
                    << rootGR << "," << QCcm << "," << avg_pressure << "," << avg_sigma << ","
                    << avg_cell_growth << "," << avg_velocity << "," << avg_auxin << "," << std_auxin << ","
                    << avg_auxinInter << "," << avg_Pin1Cyt << "," << avg_Pin1Mem << "," << avg_Aux1Mem << ","
                    << chemicalsProcess->debugs["Average Auxin Diffusion"] << ","
                    << chemicalsProcess->debugs["Average Auxin Export"]*-1 << ","
                    << chemicalsProcess->debugs["Average PIN1 Expression"] << ","
                    << chemicalsProcess->debugs["Average PIN1 Trafficking"] << ","
                    << anisotropy_degree
                    << endl
                    << flush;


        // positions, auxin and growth rate (for plotting)                
        /*
        for(auto c : cellAttr) {
            Tissue::CellData& cD = cellAttr[c.first];
            if(cD.type != Tissue::Source && cD.type != Tissue::Substrate && cD.type != Tissue::QC )
                cerr <<  mechanicsProcess->userTime << "," << cD.type << "," << cD.centroid.y() - VIcm.y() << "," << cD.auxin/cD.area << "," << cD.growthRate << "," << norm(cD.a1) << "," <<  norm(cD.a2) << endl;
        }
        */

    }

    // Calculate FPS
    clock_t curr_clock = clock();
    double elapsed_secs = double(curr_clock - prev_clock) / CLOCKS_PER_SEC;
    if(stepCount % 10 == 0)
        mdxInfo << "FPS: " << (stepCount - prevStepCount) / elapsed_secs << ", time in seconds: " <<  double(curr_clock - begin_clock) / CLOCKS_PER_SEC <<  ", steps: " << stepCount << " steps" << endl;
    prev_clock = curr_clock;
    prevStepCount = stepCount;

    // Check if execution is finished
    if(double(curr_clock - begin_clock) / CLOCKS_PER_SEC > parm("Execution Time").toDouble() && parm("Execution Time").toDouble() > 0) {
            mdxInfo << "The execution reached the maximum amount of seconds set" << endl;
            saveMeshProcess->run(mesh, parm("Output Mesh"), false);
            exit(0);
        }

    return true;
}


bool ClearCells::clearCell_old(int label) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("Root::ClearCells::clearCell No current mesh"));
    if(!getProcess(parm("AddFace"), addFaceProcess))
        throw(QString("Root::ClearCells Cannot make AddFace") + parm("AddFace"));

    CCStructure& cs = mesh->ccStructure("Tissue");
    CCIndexDataAttr& indexAttr = mesh->indexAttr();

    std::set<CCIndex> to_delete_faces;
    std::set<CCIndex> to_delete_edges;
    std::set<CCIndex> to_delete_vertices;

    std::vector<CCIndex> faces;

    updateGeometry(cs, indexAttr);

    // get faces from that label
    Point3d centroid;
    for(CCIndex f : cs.faces())
        if(indexAttr[f].label == label) {
            faces.push_back(f);
            centroid += indexAttr[f].pos;
    }
    centroid /= faces.size();

    // cell already cleared
    if(faces.size() < 2)
        return false;

    // set of edges, unorderded
    std::set<CCIndex> bn;
    for(CCIndex f : faces) {
        // ignore non-triangular faces
        /*if(cs.incidentCells(f, 1).size() != 3)
            throw(QString("Root::ClearCells No triangular cell found " + cD.label));*/
        std::set<int> labels;
        for(CCIndex e : cs.incidentCells(f, 1)) {
            labels.clear();
            for(CCIndex fn : cs.incidentCells(e, 2))
                labels.insert(indexAttr[fn].label);
            if (labels.size() > 1 || cs.onBorder(e)) { // it's a perimeter edge
                bn.insert(e);
            } else if (labels.size() == 1)
                to_delete_edges.insert(e);
        }
        for(CCIndex v : cs.incidentCells(f, 0)) {
            labels.clear();
            for(CCIndex fn : cs.incidentCells(v, 2))
                labels.insert(indexAttr[fn].label);
            if (labels.size() > 1 || cs.onBorder(v)) { // it's a perimeter vertex
               ;
            } else
                to_delete_vertices.insert(v);
        }
        to_delete_faces.insert(f);
    }
    if(bn.size() < 3)
        throw(QString("Root::ClearCells Error, not enough perimeter edges found for cell " + label));

    /*

    // order the vertices
    std::vector<std::pair<Point3d, Point3d>> polygonSegs;
    std::vector<CCIndex> vs_orig;
    std::vector<CCIndex> vs_ordered;

    // get the original edges
    for(auto i : bn) {
        std::pair<CCIndex, CCIndex> eb = cs.edgeBounds(i);
        vs_orig.push_back(eb.first);
        vs_orig.push_back(eb.second);
        polygonSegs.push_back(make_pair(indexAttr[eb.first].pos, indexAttr[eb.second].pos));
    }

    // sort them
    for(auto i : orderPolygonSegs(polygonSegs))
        for(auto j : vs_orig)
            if(i == indexAttr[j].pos) {
                vs_ordered.push_back(j);
                break;
            }

    // vertices must be ordered counter-clockwise for the face to be displayed
    // in front of the viewer
    Point3d cross_product = (indexAttr[vs_ordered[0]].pos - centroid)
                                .cross(indexAttr[vs_ordered[1]].pos - centroid);
    if(cross_product.z() < 0)
        std::reverse(std::begin(vs_ordered), std::end(vs_ordered));
    face_to_create = vs_ordered;

    */

    // delete marked edges, faces and vertices
    for(auto i : to_delete_faces)
        cs.deleteCell(i);
    for(auto i : to_delete_edges)
        cs.deleteCell(i);
    for(CCIndex v : cs.vertices())
        if(cs.cobounds(v).size() == 0)
            to_delete_vertices.insert(v);
    for(auto i : to_delete_vertices)
        cs.deleteCell(i);

    // get vertices to be used to create the new empty face
    std::set<CCIndex> vs;
    for(CCIndex e : bn) {
        vs.insert(cs.edgeBounds(e).first);
        vs.insert(cs.edgeBounds(e).second);
    }

    addFaceProcess->addFace(cs, indexAttr, vs, label);

    /*
    // add final face
    CCIndex ft = CCIndexFactory.getIndex();
    addFace(cs, ft, face_to_create);
    indexAttr[ft].label = label;
    */


    // is this necessary?
    mesh->updateAll("Tissue");

    return true;
}

CCIndex ClearCells::clearCell(int label) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("CellTissueProcess::initialize No current mesh"));

    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    CCStructure& cs =mesh->ccStructure("Tissue");
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");

    CCIndexTbbSet delFaces, delEdges, delVertices;

    for(CCIndex f: *(cellAttr[label].cellFaces))
        delFaces.insert(f);

    // Leave if no cells selected for deletion
    if(delFaces.size() == 0)
      return CCIndex();

    // Delete the dangling edges
    for(CCIndex f : delFaces) {
        for(CCIndex e : cs.incidentCells(f,1))
            if(find(cellAttr[label].perimeterEdges.begin(), cellAttr[label].perimeterEdges.end(), e) == cellAttr[label].perimeterEdges.end())
                  delEdges.insert(e);
    }

    // Delete the dangling vertices
    for(CCIndex f : delFaces) {
        for(CCIndex v : cs.incidentCells(f,0))
            if(find(cellAttr[label].perimeterVertices.begin(), cellAttr[label].perimeterVertices.end(), v) == cellAttr[label].perimeterVertices.end())
                  delVertices.insert(v);
    }

    for(CCIndex f : delFaces)
        cs.deleteCell(f);

    for(CCIndex e : delEdges)
        cs.deleteCell(e);

    for(CCIndex v : delVertices)
        cs.deleteCell(v);

    CCIndex ft = CCIndexFactory.getIndex();
    if(!mdx::addFace(cs, ft, cellAttr[label].perimeterVertices)) {
        std::reverse(cellAttr[label].perimeterVertices.begin(), cellAttr[label].perimeterVertices.end());
        mdx::addFace(cs, ft, cellAttr[label].perimeterVertices);
    }
    indexAttr[ft].label = label;

    mesh->updateAll();

    return ft;

}


bool ClearCells::clearAllCells() {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("Root::ClearCells No current mesh"));

    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");

    mdxInfo << "ClearCells::clearAllCells: Clearing " << cellAttr.size() << " cells" << endl;
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        clearCell(cD.label);
    }

    return false;
}


bool TriangulateFacesX::step_destructive(CCStructure& cs,
                             CCStructure& csOut,
                             CCIndexDataAttr& indexAttr,
                             double maxArea) {

    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("TriangulateFaces No current mesh"));


    mdxInfo << "TriangulateFacesX::step_destructive" << endl;

    Splitter subdiv;
    // Setup subdivision objects
    subdiv.mdx = MDXSubdivideX(*mesh);
    subdiv.mechanics =
        Tissue::Subdivide(indexAttr,
                          mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData"),
                          mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData"),
                          mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData"),
                          mesh->attributes().attrMap<int, Tissue::CellData>("CellData"));


    // First split edges
    double maxLength = sqrt(2 * maxArea);
    do {
        CCIndexVec edges;
        for(CCIndex e : cs.edges()) {
            auto eb = cs.edgeBounds(e);
            if(norm(indexAttr[eb.first].pos - indexAttr[eb.second].pos) > maxLength)
                edges.push_back(e);
        }
        if(edges.size() > 0)
            splitEdges(cs, indexAttr, edges, &subdiv);
        else
            break;
    } while(true);

    bool boundary = false;
    bool inOrder = true;
    bool center = true;

    int dim = cs.maxDimension;
    MeshBuilder mb(indexAttr, cs.maxDimension);

    // triangulate using triangle (jajajajaj thank you captain obvious)
    for(CCIndex f : cs.faces()) {
        auto& fIdx = indexAttr[f];

        CCIndexVec fVertices = faceVertices(cs, f);
        Point3dVec poly(fVertices.size());
        for(uint i = 0; i < fVertices.size(); i++)
            poly[i] = indexAttr[fVertices[i]].pos;

        Point3d polyNormal;
        std::vector<Point3i> triList;
        Point3dVec ptList;

        // calculate the polygon plane using pca
        Point3d planePos, planeNrml;
        Point3dVec polyProj = poly;
        findPolygonPlane(polyProj, planePos, planeNrml);
        // Why does this seem backward?
        if(fIdx.nrml * planeNrml > 0)
            planeNrml = -planeNrml;

        // project all points onto the plane
        projectPointsOnPlane(polyProj, planePos, planeNrml);

        // triangulate the polygon
        if(!triangulatePolygon3D(
               maxArea, planeNrml, polyProj, triList, ptList, boundary, inOrder, center))
            mdxInfo << QString("%1: Error calling triangulatePolygon3D").arg(name()) << endl;

        // move the border points back to their original position
        if(ptList.size() < poly.size()) {
            mdxInfo << QString("%1: Error in ptList: %2 poly: %3 triList %4")
                           .arg(name())
                           .arg(ptList.size())
                           .arg(triList.size())
                           .arg(poly.size())
                    << endl;
            continue;
        }
        for(uint i = 0; i < poly.size(); i++)
            ptList[i] = poly[i];

        // create the triangular faces
        if(dim == 2) {
            int label = indexAttr[f].label;
            for(const Point3i& tri : triList)
                mb.addFace({ptList[tri[0]], ptList[tri[1]], ptList[tri[2]]}, label);
        } else if(dim == 3) {
            auto cobounds = cs.cobounds(f);
            int label = 0, nlabel = 0;
            if(cobounds.size() > 0)
                label = indexAttr[*cobounds.begin()].label;
            if(cobounds.size() > 1)
                nlabel = indexAttr[*cobounds.rbegin()].label;
            if(cs.ro(*cobounds.begin(), f) == ccf::NEG) {
                if(nlabel == 0)
                    mdxInfo << QString("%1: Warning, outside face has the wrong orientation")
                                   .arg(name())
                            << endl;
                else
                    std::swap(label, nlabel);
            }
            if(label > 0)
                for(const Point3i& tri : triList)
                    mb.addFace({ptList[tri[0]], ptList[tri[1]], ptList[tri[2]]}, label, nlabel);
            else if(nlabel > 0)
                for(const Point3i& tri : triList)
                    mb.addFace({ptList[tri[0]], ptList[tri[2]], ptList[tri[1]]}, nlabel);
        }
    }

    csOut = CCStructure(dim);
    mb.createMesh(csOut);
    mesh->updateAll();

    return true;
}


bool TriangulateFacesX::triangulateFace(CCIndex f, CCStructure& cs,
                                        CCIndexDataAttr& indexAttr,
                                        double maxArea) {

    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("TriangulateFaces No current mesh"));


    Splitter subdiv;
    // Setup subdivision objects
    subdiv.mdx = MDXSubdivideX(*mesh);
    subdiv.mechanics =
        Tissue::Subdivide(indexAttr,
                          mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData"),
                          mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData"),
                          mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData"),
                          mesh->attributes().attrMap<int, Tissue::CellData>("CellData"));

    if(f.isPseudocell())
        throw(QString("TriangulateFacesX::triangulateFace: no face of label " + indexAttr[f].label));

    // First split edges
    std::set<CCIndex> splitted_edges;
    double maxLength = sqrt(2 * maxArea);
    do {
        std::vector<CCIndex> edges;
        for(CCIndex e : cs.edges()) {
            if(!cs.incident(f, e))
                continue;
            auto eb = cs.edgeBounds(e);
            if(norm(indexAttr[eb.first].pos - indexAttr[eb.second].pos) > maxLength) {
                splitted_edges.insert(e);
                edges.push_back(e);
            }
        }
        if(edges.size() > 0)
            splitEdges(cs, indexAttr, edges, &subdiv);
        else
            break;
    } while(true);

    bool boundary = false;
    bool inOrder = true;
    bool center = true;

    int dim = cs.maxDimension;
    std::map<Point3d, CCIndex> vertices_map;
    for(CCIndex v : cs.vertices())
       vertices_map[indexAttr[v].pos] = v;

    auto& fIdx = indexAttr[f];
    CCIndexVec fVertices = faceVertices(cs, f);
    /*if(fVertices.size() == 3)
        mdxInfo << "Nothing to do for face " << indexAttr[f].label << endl;*/
    Point3dVec poly(fVertices.size());
    for(uint i = 0; i < fVertices.size(); i++)
        poly[i] = indexAttr[fVertices[i]].pos;

    std::vector<Point3i> triList;
    Point3dVec ptList;

    // calculate the polygon plane using pca
    Point3d planePos, planeNrml;
    Point3dVec polyProj = poly;
    findPolygonPlane(polyProj, planePos, planeNrml);
    // Why does this seem backward?
    if(fIdx.nrml * planeNrml > 0)
        planeNrml = -planeNrml;

    // project all points onto the plane
    projectPointsOnPlane(polyProj, planePos, planeNrml);

    // triangulate the polygon
    if(!triangulatePolygon3D(
           maxArea, planeNrml, polyProj, triList, ptList, boundary, inOrder, center))
        mdxInfo << QString("%1: Error calling triangulatePolygon3D").arg(name()) << endl;

    // move the border points back to their original position
    if(ptList.size() < poly.size()) {
        mdxInfo << QString("%1: Error in ptList: %2 poly: %3 triList %4")
                       .arg(name())
                       .arg(ptList.size())
                       .arg(triList.size())
                       .arg(poly.size())
                << endl;
        return false;
    }
    for(uint i = 0; i < poly.size(); i++)
        ptList[i] = poly[i];

    // create the triangular faces
    if(dim == 2) {
        int label = indexAttr[f].label;
        cs.deleteCell(f);
        for(const Point3i& tri : triList) {
            CCIndex fi = CCIndexFactory.getIndex();
            indexAttr[fi].label = label;
            CCIndex v1;
            if(vertices_map.find(ptList[tri[0]]) == vertices_map.end()) {
                v1 = CCIndexFactory.getIndex();
                vertices_map[ptList[tri[0]]] = v1;
            }
            else
                v1 = vertices_map[ptList[tri[0]]];
            CCIndex v2;
            if(vertices_map.find(ptList[tri[1]]) == vertices_map.end()) {
                v2 = CCIndexFactory.getIndex();
                vertices_map[ptList[tri[1]]] = v2;
            }
            else
                v2 = vertices_map[ptList[tri[1]]];
            CCIndex v3;
            if(vertices_map.find(ptList[tri[2]]) == vertices_map.end()) {
                v3 = CCIndexFactory.getIndex();
                vertices_map[ptList[tri[2]]] = v3;
            }
            else
                v3 = vertices_map[ptList[tri[2]]];
            std::vector<CCIndex> vs;
            vs.push_back(v1);
            vs.push_back(v2);
            vs.push_back(v3);
            addFace(cs, fi, vs);
            indexAttr[v1].pos = ptList[tri[0]];
            indexAttr[v2].pos = ptList[tri[1]];
            indexAttr[v3].pos = ptList[tri[2]];
            indexAttr[fi].label = label;
            //mb.addFace({ptList[tri[0]], ptList[tri[1]], ptList[tri[2]]}, label);
        }
    } else if(dim == 3) {
        auto cobounds = cs.cobounds(f);
        int label = 0, nlabel = 0;
        if(cobounds.size() > 0)
            label = indexAttr[*cobounds.begin()].label;
        if(cobounds.size() > 1)
            nlabel = indexAttr[*cobounds.rbegin()].label;
        if(cs.ro(*cobounds.begin(), f) == ccf::NEG) {
            if(nlabel == 0)
                mdxInfo << QString("%1: Warning, outside face has the wrong orientation")
                               .arg(name())
                        << endl;
            else
                std::swap(label, nlabel);
        }
        if(label > 0) {}
            /*for(const Point3i& tri : triList){
                //////////////////////////////////////////////////////// UNSET copy from above in dim == 2
                //mb.addFace({ptList[tri[0]], ptList[tri[1]], ptList[tri[2]]}, label, nlabel);
            }*/
        else if(nlabel > 0) {}
            /*for(const Point3i& tri : triList){
                //////////////////////////////////////////////////////// UNSET, copy from above in dim == 2
                //mb.addFace({ptList[tri[0]], ptList[tri[2]], ptList[tri[1]]}, nlabel);
            }*/
    }

    mesh->updateAll();

    return true;
}

bool TriangulateFacesX::step_nondestructive(CCStructure& cs,
                             CCIndexDataAttr& indexAttr,
                             double maxArea) {

    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("TriangulateFaces No current mesh"));

    int round = 0;
    int triangulated = 0;
    do {
        triangulated = 0;
        CCIndexVec fs = cs.faces();
        for(CCIndex f : fs) {
            if(cs.incidentCells(f, 1).size() > 3) {
                triangulateFace(f, cs, indexAttr, maxArea);
                triangulated++;
            }
        }
        round ++;
        //mdxInfo << "TriangulateFacesX::step_nondestructive: I have triangulated " <<  triangulated << " faces in round " << round << endl;
    } while(triangulated);

    mesh->updateAll();

    return true;
}


bool PrintCellAttr::step() {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("PrintCellAttr::step No current mesh, cannot rewind"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
        throw(QString("PrintCellAttr::run Error, no cell complex selected"));

    CCStructure& cs = mesh->ccStructure("Tissue");
    CCStructure& csDual = mesh->ccStructure("TissueDual");

    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    Tissue::FaceDataAttr& faceAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");
    mdxInfo << endl << endl;
    for(CCIndex f : cs.faces()) {
        CCIndexData& cIdx = indexAttr[f];
        int label = indexAttr[f].label;
        if(cIdx.selected) {
            Tissue::CellData cD = cellAttr[label];
            mdxInfo << "Cell owner: " << cD.label << endl
                    << " type: " << Tissue::ToString(cD.type) << endl;
            mdxInfo << " perimeter faces: ";
            for(auto f : cD.perimeterFaces)
                mdxInfo << f << " ";
            mdxInfo << endl;
            mdxInfo << " perimeter edges: ";
            for(auto e : cD.perimeterEdges)
               mdxInfo << e << " ";// << " MF Impact " << edgeAttr[e].MFImpact[label] << " auxinflux " << edgeAttr[e].auxinFluxImpact[label] << " raw " <<  edgeAttr[e].pin1SensitivityRaw[label] << " norm " << edgeAttr[e].pin1Sensitivity[label] << endl ;

            mdxInfo << endl;
            mdxInfo << " perimeter vertices: ";
            for(auto v : cD.perimeterVertices)
                mdxInfo << v << " angle: " << cD.perimeterAngles[v]*(180/M_PI) << " ";
            mdxInfo << endl;
            mdxInfo << " my neighbours are: ";
            if(csDual.hasCell(cD.dualVertex))
                for(CCIndex vn : csDual.neighbors(cD.dualVertex)) {
                    int labeln = indexAttr[vn].label;
                    Tissue::CellData& cDn = cellAttr[labeln];
                    mdxInfo << cDn.label << ", " ;
                }
            mdxInfo << endl;
            mdxInfo
                    << " center: " << cD.centroid << " dualVertex " << cD.dualVertex <<  endl
                    << " area: " << cD.area << " "
                    << " restArea: " << cD.restArea << " "
                    << " prevArea: " << cD.prevArea << " "
                    << " max area: " << cD.cellMaxArea << " "
                    << " perimeter: " << cD.perimeter << " "
                    << " invmassVertices: " << cD.invmassVertices << " " << endl
                    << " G: " << cD.G << " E: " << cD.E << " F: " << cD.F << " R: " << cD.R << " U: " << cD.U << " S: " << cD.S << " M: " << cD.M << " M0: " << cD.M0
                    << " gMax: " << cD.gMax
                    << " gMin: " << cD.gMin
                    << " restCm: " << cD.restCm << " "
                    << " invRestMat: " << cD.invRestMat << " "
                    << " restX0: ";
                    for(auto x : cD.restX0) mdxInfo << x << " ";
                    mdxInfo << endl;
            mdxInfo << " a1: " << cD.a1 << " " << " a2: " << cD.a2 << " "
                    << " axisMin: " << cD.axisMin << " " << " axisMax: " << cD.axisMax << " " << " divVector " << cD.divVector << " MF reorientation: " << cD.mfRORate << endl
                    << " periclinal division: " << cD.periclinalDivision <<  " division algorithm: " << cD.divAlg << " last division: " << cD.lastDivision << " division counts: " << cD.divisionCount << endl
                    << " pressure: " << cD.pressure << " " << " pressureMax: " << cD.pressureMax << " " << " wallStress : " << cD.wallStress << endl;
            mdxInfo << " strain rate on the edges: ";
                       for(auto e : cD.perimeterEdges) {
                           Tissue::EdgeData& eD = edgeAttr[e];
                           mdxInfo << eD.strainRate << " ";
                       }
            mdxInfo << endl;
            mdxInfo
                    << " auxin: " << cD.auxin << " " << " auxin by area: " << cD.auxin/cD.area << " "
                    //<< " growth factor: " << cD.growthFactor << " " << " growth factor by area: " << cD.growthFactor/cD.area << " "
                    << " Aux1: " << cD.Aux1 << " "
                    << " Pin1: " << cD.Pin1 << " " << " Pin1 by area: " << cD.Pin1/cD.area << " "
                    << " Division Promoter: " << cD.divPromoter/cD.area << " " << " Division Inhibitor: " << cD.divInhibitor/cD.area << " "<< " Division Probability: " << cD.divProb << " "
                    << " PINOID: " << cD.PINOID << " "   << " PP2A: " << cD.PP2A << " "
                    << " pinProdRate: " << cD.pinProdRate << " " << " aux1ProdRate: " << cD.aux1ProdRate << " "<< " pinInducedRate: " << cD.pinInducedRate << " " << " aux1InducedRate: " << cD.aux1InducedRate << " "<< " aux1MaxEdge: " << cD.aux1MaxEdge << " "
                    << " auxinProdRate: " << cD.auxinProdRate << " " << " auxinFluxVector: " << cD.auxinFluxVector
                    << endl;
            for(auto p : cD.auxinFluxes)
                mdxInfo << " auxinFlux with " << p.first << " : " << p.second << " ";
            mdxInfo << endl << endl;
            Tissue::FaceData fD = faceAttr[f];
            mdxInfo << "Face: " << f  << " type: " << Tissue::ToString(fD.type) << " center: " << cIdx.pos << " label " << cIdx.label  << " owner " << fD.owner
                    << endl
                    //<< " orientation with TOP: " << cs.ro(CCIndex::TOP, f) << endl
                    << " area: " << fD.area << " restAreaFace: " << fD.restAreaFace << " invRestMat: " << fD.invRestMat << " restPos: " << fD.restPos[0] << " : " << fD.restPos[1] << " : " << fD.restPos[2]
                    << endl
                    //<< " prev_centroid: " << fD.prev_centroid << endl
                    << " E: " << fD.E << " " << " F: " << fD.F << " " << " R: " << fD.R << " " << " G: " << fD.G << " " << " V: " << fD.V << endl
                    //<< " F1: " << fD.F1 << " " << " F2: " << fD.F2 << endl
                    << " a1: " << fD.a1 << " " << " a2: " << fD.a2 << endl
                    << " stress: " << fD.stress << " " << " sigmaA: " << fD.sigmaA << endl
                    << " (visuals)"
                    << " type : " << fD.type << " growthRate : " << fD.growthRate << " auxin : " << fD.auxin << " intercellularAuxin : " << fD.intercellularAuxin
                    << " Pin1Cyt : " << fD.Pin1Cyt << " Pin1Mem : " << fD.Pin1Mem
                    << " Aux1Cyt : " << fD.Aux1Cyt << " Aux1Mem : " << fD.Aux1Mem
                    << endl
                    << endl;
            for(auto e : cs.incidentCells(f, 1)) {
                Tissue::EdgeData& eD = edgeAttr[e];
                mdxInfo << "Edge: " << e << " " << indexAttr[cs.edgeBounds(e).first].pos << " : "
                        << indexAttr[cs.edgeBounds(e).second].pos << " type: " << Tissue::ToString(eD.type)
                        << " length: " << eD.length << " prevLength: " << eD.prevLength
                        << " restLength: " << eD.restLength ;
                for(auto p : eD.outwardNormal)
                   mdxInfo << " outward normal to " << indexAttr[p.first].label << " : " << p.second;
                mdxInfo << endl ;
                mdxInfo << " cStiffness: " << eD.cStiffness << " eStiffness: " << eD.eStiffness
                        << " prev_strain: " << eD.prev_strain
                        << " strain: " << eD.strain << " strainRate: " << eD.strainRate
                        //<< " cellAxis: " << eD.cellAxis
                        << " sigmaEv: " << eD.sigmaEv << " sigmaEe: " << eD.sigmaEe << " sigmaForce: " << eD.sigmaForce
                        << " pressureForce: " << eD.pressureForce  << endl;
                 mdxInfo  << " intercellularAuxin: " << eD.intercellularAuxin << endl;
                 for(auto p : eD.auxinRatio)
                     mdxInfo << " auxinRatio from " << p.first << " : " << p.second;
                 mdxInfo << endl ;
                 for(auto p : eD.auxinGrad)
                     mdxInfo << " auxinGrad from " << p.first << " : " << p.second;
                 mdxInfo << endl ;
                for(auto p : eD.MFImpact)
                    mdxInfo << " MFImpact from " << p.first << " : " << p.second;
                mdxInfo << endl ;
                for(auto p : eD.auxinFluxImpact)
                    mdxInfo << " auxinFluxImpact from " << p.first << " : " << p.second;
                mdxInfo << endl ;
                for(auto p : eD.geomImpact)
                    mdxInfo << " geommpact from " << p.first << " : " << p.second;
                mdxInfo << endl ;
                for(auto p : eD.pin1SensitivityRaw)
                    mdxInfo << " pin1SensitivityRaw from " << p.first << " : " << p.second;
                for(auto p : eD.pin1Sensitivity)
                    mdxInfo << " pin1Sensitivity from " << p.first << " : " << p.second;
                mdxInfo << endl ;
                for(auto p : eD.Pin1)
                    mdxInfo << " Pin1 from " << p.first << " : " << p.second;
                mdxInfo << endl;
                for(auto p : eD.Aux1)
                    mdxInfo << " Aux1 from " << p.first << " : " << p.second;
                mdxInfo << endl;
                for(auto p : eD.PINOID)
                    mdxInfo << " PINOID from " << p.first << " : " << p.second;
                mdxInfo << endl;
                for(auto p : eD.PP2A)
                    mdxInfo << " PP2A from " << p.first << " : " << p.second;
                mdxInfo << endl;
            }
            mdxInfo << endl;
            for(auto v : cs.incidentCells(f, 0)) {
                Tissue::VertexData vD = vMAttr[v];
                mdxInfo << "Vertex:" << v << " " << indexAttr[v].pos
                        << " type: " << Tissue::ToString(vD.type)
                        << " invmass: " << vD.invmass;
                        /*for(auto a : vD.angle)
                            mdxInfo << " angle with " << a.first[f] << ": " << a.second * (180/M_PI);*/
                        mdxInfo << " dualCell: " << vD.dualCell
                        << " divPoint: " << vD.divisionPoint
                        << " prevPos: " << vD.prevPos
                        << " restPos: " << vD.restPos
                        << " velocity: " << vD.velocity
                        << " prevVelocity: " << vD.prevVelocity
                        << " norm velocity: " << norm(vD.velocity)
                        << " dampedVelocity: " << norm(vD.dampedVelocity)
                        << " forces: " << vD.forces.size();
                Point3d pressure, sigmaAY, sigmaEv, sigmaEe, totalForce;
                for(auto m : vD.forces) {
                    mdxInfo << " , " << std::get<0>(m) << " " << std::get<1>(m) << " "
                            << std::get<2>(m);
                    totalForce += std::get<2>(m);
                    if(std::get<1>(m) == QString("pressure"))
                        pressure += std::get<2>(m);
                    else if(std::get<1>(m) == QString("sigmaAY"))
                        sigmaAY += std::get<2>(m);
                    else if(std::get<1>(m) == QString("sigmaEv"))
                        sigmaEv += std::get<2>(m);
                    else if(std::get<1>(m) == QString("sigmaEe"))
                        sigmaEe += std::get<2>(m);
                    else if(std::get<1>(m) == QString("substrate") ||
                            std::get<1>(m) == QString("damping") ||
                            std::get<1>(m) == QString("friction"))
                        ;
                    else
                        throw(QString("PrintCellAttr::step Unknown force: " +
                                      QString(std::get<1>(m))));
                }
                mdxInfo << endl;
                mdxInfo << " Total Force: " << totalForce << endl
                        << "  Pressure: " << pressure
                        << "  sigmaAY: " << sigmaAY
                        << "  sigmaEv: " << sigmaEv
                        << "  sigmaEe: " << sigmaEe << endl;;
                Point3d totalCorrection;
                for(auto c : vD.corrections)
                    totalCorrection += c.second;
                mdxInfo << " Total Correction: " << totalCorrection << endl;
                for(auto c : vD.corrections)
                    mdxInfo << "  " << c.first << ": " << c.second;
                mdxInfo << endl;

            }
        }
    }
    mdxInfo << "Total cells: " << cellAttr.size() << endl;
    /*mdxInfo << "Cell set: " << endl;
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        mdxInfo << "Label: " << cD.label << " pos: " << cD.centroid << endl;
    }*/
    return false;
}


bool SetGlobalAttr::initialize(QWidget* parent) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("SetGlobalAttr::initialize No current mesh"));
    if(!getProcess(parm("Mechanics Process"), mechanicsProcess))
        throw(QString("SetGlobalAttr::initialize Cannot make Mechanics Process:") + parm("Mechanics Process"));
    if(!getProcess(parm("Divide Process"), divideProcess))
        throw(QString("SetGlobalAttr::initialize Cannot make Divide Process:") + parm("Divide Process"));
    return step();
}

bool SetGlobalAttr::step() {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("SetCell::step No current mesh"));

    CCStructure& cs = mesh->ccStructure("Tissue");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");

    mdxInfo << "Setting global attrs" << endl;

    // wall edges are shared between cells so we will take the average
    for(CCIndex e : cs.edges())
        if(edgeAttr[e].type == Tissue::Wall)
            edgeAttr[e].cStiffness = edgeAttr[e].eStiffness = 0;
    for(auto c : cellAttr) {  // If values != -1 all paramters are replaced
        Tissue::CellData& cD = cellAttr[c.first];
        QString str = Tissue::ToString(cD.type);
        if(parm(QString(str + " Turgor Pressure")).toDouble() >= 0)
            cD.pressureMax = parm(QString(str + " Turgor Pressure")).toDouble();
        cD.growthFactor = parm(QString(str + " Growth Factor")).toDouble();
        if(parm(QString(str + " Max area")).toDouble() >= 0)
            cD.cellMaxArea = parm(QString(str + " Max area")).toDouble();
        else
            cD.cellMaxArea = divideProcess->parm("Cell Max Area").toDouble();
        if(parm(QString(str + " MF reorientation rate")).toDouble() > -1)
            cD.mfRORate = parm(QString(str + " MF reorientation rate")).toDouble();
        for(CCIndex f : *cD.cellFaces)
            for(CCIndex e : cs.incidentCells(f, 1)) {
                Tissue::EdgeData& eD = edgeAttr[e];
                double wallCK = parm(QString(str + " Wall CK")).toDouble();
                double wallEK = parm(QString(str + " Wall EK")).toDouble();
                double shearCK = parm(QString(str + " Shear CK")).toDouble();
                double shearEK = parm(QString(str + " Shear EK")).toDouble();
                if(eD.type == Tissue::Wall) {
                    if(wallCK >= 0)
                        eD.cStiffness += wallCK;
                    else
                        eD.cStiffness += mechanicsProcess->parm("Wall CK").toDouble();
                    if(wallEK >= 0)
                        eD.eStiffness += wallEK;
                    else
                        eD.eStiffness += mechanicsProcess->parm("Wall EK").toDouble();

                }
                else if(eD.type == Tissue::Shear) {
                    if(shearCK >= 0)
                        eD.cStiffness = shearCK;
                    else
                        eD.cStiffness = mechanicsProcess->parm("Shear CK").toDouble();
                    if(shearEK >= 0)
                        eD.eStiffness = shearEK;
                    else
                        eD.eStiffness = mechanicsProcess->parm("Shear EK").toDouble();
                }
            }
    }
    for(CCIndex e : cs.edges())
        if(edgeAttr[e].type == Tissue::Wall)
            if(!cs.onBorder(e)) {
                edgeAttr[e].cStiffness /= 2;
                edgeAttr[e].eStiffness /= 2;
            }
    return false;
}

bool SetCellAttr::step() {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("SetGlobalAttr::step No current mesh, cannot rewind"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
        throw(QString("SetGlobalAttr::step Error, no cell complex selected"));

    CCStructure& cs = mesh->ccStructure(ccName);
    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    Tissue::FaceDataAttr& faceAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");

    for(CCIndex f : cs.faces()) {
        CCIndexData& cIdx = indexAttr[f];
        if(cIdx.selected) {
            Tissue::CellData& cD = cellAttr[cIdx.label];
            cD.type = Tissue::stringToCellType(parm("Cell Type"));
            cD.divAlg = parm("Division Algorithm").toInt();
            cD.auxinProdRate = parm("Auxin production rate").toDouble();
            if(parm("Periclinal Division") != "")
                cD.periclinalDivision = parm("Periclinal Division") == "True";
            cD.mfRORate = parm("MF reorientation rate").toDouble();
            QStringList list = parm("Microfibril 1").split(QRegExp(","));
            cD.a1[0] = list[0].toDouble();
            cD.a1[1] = list[1].toDouble();
            cD.a1[2] = list[2].toDouble();
            list = parm("Microfibril 2").split(QRegExp(","));
            cD.a2[0] = list[0].toDouble();
            cD.a2[1] = list[1].toDouble();
            cD.a2[2] = list[2].toDouble();
            double mass =  parm("Vertices Masses").toDouble();
            cD.invmassVertices = 1. / mass;
            for(CCIndex f : (*cD.cellFaces)) {
                Tissue::FaceData& fD = faceAttr[f];
                fD.a1 = cD.a1;
                fD.a2 = cD.a2;
                for(CCIndex v : cs.incidentCells(f, 0))
                    vMAttr[v].invmass = cD.invmassVertices;
            }
            cD.restX0.clear();
            std::vector<double> invMasses;
            for(CCIndex v : cD.cellVertices) {
                cD.restX0.push_back(Point3d(indexAttr[v].pos));
                invMasses.push_back(vMAttr[v].invmass);
            }
            PBD::init_ShapeMatchingConstraint(cD.restX0, invMasses, cD.restX0.size(),  cD.restCm, cD.invRestMat);

        }
    }

    return false;
}

bool HighlightCell::step() {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("HighlightCell::step No current mesh"));
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");

    CCStructure& cs = mesh->ccStructure("Tissue");
    CCIndexDataAttr& indexAttr = mesh->indexAttr();

    for(CCIndex f : cs.faces())
        if(indexAttr[f].label == parm("Face Label").toInt())
            indexAttr[f].selected = true;
    if(cellAttr.find(parm("Face Label").toInt()) == cellAttr.end())
        throw(QString("Cell does not exists"));
    mesh->updateAll("Tissue");


    return false;
}





REGISTER_PROCESS(Root);
REGISTER_PROCESS(Tissue);
REGISTER_PROCESS(PBD);
REGISTER_PROCESS(Mechanics);
REGISTER_PROCESS(MechanicalGrowth);
REGISTER_PROCESS(Chemicals);
REGISTER_PROCESS(Remesh);
REGISTER_PROCESS(RootDivide);
REGISTER_PROCESS(SetCellAttr);
REGISTER_PROCESS(SetGlobalAttr);
REGISTER_PROCESS(PrintCellAttr);
REGISTER_PROCESS(HighlightCell);
REGISTER_PROCESS(DeleteCell);
REGISTER_PROCESS(ExecuteTest);
REGISTER_PROCESS(ClearCells);
REGISTER_PROCESS(TriangulateFacesX);
REGISTER_PROCESS(PrintVertexAttr);
// REGISTER_PROCESS(SplitEdges);
REGISTER_PROCESS(AddFace);
REGISTER_PROCESS(DeleteEdges);
REGISTER_PROCESS(CreateEdge);
REGISTER_PROCESS(ReverseCell);
