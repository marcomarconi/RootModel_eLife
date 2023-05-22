
#include "tissue.hpp"


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%T", &tstruct);

    return buf;
}

void SVDDecompX(const Matrix3d &M, Matrix3d &U, Matrix3d &S, Matrix3d &V){
   gsl_matrix* Mgsl = gsl_matrix_alloc(3, 3);
   gsl_vector* Sgsl = gsl_vector_alloc(3);
   gsl_matrix* Vgsl = gsl_matrix_alloc(3, 3);
   gsl_vector* wgsl = gsl_vector_alloc(3);


   for(int i = 0; i < 3; ++i)
     for(int j = 0; j < 3; ++j)
       gsl_matrix_set(Mgsl, i, j, M[i][j]);
   gsl_linalg_SV_decomp(Mgsl, Vgsl, Sgsl, wgsl);

   for(int i = 0; i < 3; ++i)
     for(int j = 0; j < 3; ++j){
        U[i][j] = gsl_matrix_get(Mgsl,i,j);
        V[i][j] = gsl_matrix_get(Vgsl,i,j);
     }
   //This could be passed around as just a 3D vector, but we'll use a matrix for now
   S = Matrix<3,3,double>::Diagonal(Point3d(gsl_vector_get(Sgsl, 0), gsl_vector_get(Sgsl, 1), gsl_vector_get(Sgsl, 2)));

   gsl_matrix_free(Mgsl);
   gsl_matrix_free(Vgsl);
   gsl_vector_free(Sgsl);
   gsl_vector_free(wgsl);
}

void PolarDecompX(Matrix3d &M,Matrix3d &S, Matrix3d &U){
  Matrix3d Usvd;
  Matrix3d Ssvd;
  Matrix3d Vsvd;

  SVDDecompX(M,Usvd,Ssvd,Vsvd);

  U = Usvd*transpose(Vsvd);
  S = Vsvd*Ssvd*transpose(Vsvd);
}




Point3d Rotate(Point3d v, double angle)
{
    Point3d p = v;
    double Tx=p[0],Ty=p[1];

    double X=(Tx*cos(angle))-(Ty*sin(angle));
    double Y=(Ty*cos(angle))+(Tx*sin(angle));

    p[0] = X;
    p[1] = Y;

    return p;
}

bool PointInPolygon(Point3d point, std::vector<Point3d> points) {
    int i, j, nvert = points.size();
    bool c = false;

    for(i = 0, j = nvert - 1; i < nvert; j = i++) {
        if(((points[i].y() > point.y()) != (points[j].y() > point.y())) &&
           (point.x() < (points[j].x() - points[i].x()) * (point.y() - points[i].y()) /
                                (points[j].y() - points[i].y()) +
                            points[i].x()))
            c = !c;
    }

    return c;
}

Point3d shortenLength(Point3d A, float reductionLength) {
    Point3d B = A;
    B *= (1 - reductionLength / A.norm());
    return B;
}

float DistancePtLine(Point3d a, Point3d b, Point3d p) {
    Point3d n = b - a;
    Point3d pa = a - p;
    Point3d c = n * ((pa * n) / (n * n));
    Point3d d = pa - c;
    return sqrt((d * d));
}

Point3d findClosestLineToLine(Point3d targetLine,
                           Point3d line1, Point3d line2) {
    Point3d finalLine;
    double minAxis1;
    if(mdx::angle(targetLine, line1) < mdx::angle(targetLine, -line1))
         minAxis1 = mdx::angle(targetLine, line1);
    else
         minAxis1 = mdx::angle(targetLine, -line1);
    double minAxis2;
    if(mdx::angle(targetLine, line2) < mdx::angle(targetLine, -line2))
         minAxis2 = mdx::angle(targetLine, line2);
    else
         minAxis2 = mdx::angle(targetLine, -line2);
    if(minAxis1 < minAxis2)
        finalLine = line1;
    else
        finalLine = line2;
    return finalLine;
}

std::vector<double> softmax(std::vector<double> v) {
    double sum = 0;
    std::vector<double> r;
    for(double i : v) {
        r.push_back(exp(i));
        sum += exp(i);
    }
    for(uint i = 0; i < r.size(); i++)
        r[i] /= sum;
    return r;
}

void MBR(std::vector<Point3d> points, std::vector<Point3d>& rect) {
    std::vector<int> permutation;
    sortPolygonPoints(points, Point3d(0, 0, -1), permutation);

    std::vector<Point3d> p;
    for(auto i : permutation)
        p.push_back(points[i]);
    p.push_back(p[0]);

    std::vector<Point3d> e;
    for(uint i = 1; i < p.size(); i++)
        e.push_back(p[i] - p[i - 1]);

    std::vector<double> norms;
    for(auto i : e)
        norms.push_back(norm(i));

    std::vector<Point3d> v;
    for(uint i = 0; i < e.size(); i++)
        v.push_back(e[i] / norms[i]);

    std::vector<Point3d> w;
    for(uint i = 0; i < v.size(); i++)
        w.push_back(Point3d(-v[i][1], v[i][0], 0));

    std::vector<double> a, b, areas;
    std::vector<std::vector<double>> x, y;
    for(uint i = 0; i < v.size(); i++) {
        a.clear();
        b.clear();
        for(uint j = 0; j < p.size(); j++) {
            a.push_back(p[j] * v[i]);
            b.push_back(p[j] * w[i]);
        }

        double max_a = *max_element(a.begin(), a.end());
        double min_a = *min_element(a.begin(), a.end());
        double max_b = *max_element(b.begin(), b.end());
        double min_b = *min_element(b.begin(), b.end());
        areas.push_back((min_b - max_b) * (min_a - max_a));

        x.push_back({min_a, max_a});
        y.push_back({min_b, max_b});
    }

    int k = std::min_element(areas.begin(), areas.end()) - areas.begin();
    int qx[4] = {0, 1, 1, 0};
    int qy[4] = {0, 0, 1, 1};
    Matrix<4, 2, double> M1;
    for(int i = 0; i < 4; i++)
        M1[i] = Point2d(x[k][qx[i]], y[k][qy[i]]);
    Matrix<2, 3, double> M2;
    M2[0] = v[k];
    M2[1] = w[k];
    Matrix<4, 3, double> M = M1 * M2;
    rect.clear();
    for(int i = 0; i < 4; i++)
        rect.push_back(M[i]);
}

// extended version of the official one
void neighborhood2D(Mesh& mesh,
                    const CCStructure& cs,
                    const CCIndexDataAttr& indexAttr,
                    std::map<IntIntPair, double>& wallAreas,
                    std::map<IntIntPair, std::set<CCIndex>>& wallEdges) {

    wallAreas.clear();
    wallEdges.clear();
    //#pragma omp parallel for
    for(uint i = 0; i < cs.edges().size(); i++) {
        CCIndex e = cs.edges()[i];
        std::set<CCIndex> edgeVtxs = cs.incidentCells(e, 0);
        CCIndex v1 = *(edgeVtxs.begin());
        CCIndex v2 = *(++edgeVtxs.begin());

        double edgeLength = norm(indexAttr[v1].pos - indexAttr[v2].pos);

        std::set<CCIndex> edgeFaces = cs.incidentCells(e, 2);
        CCIndex f1 = *(edgeFaces.begin());
        CCIndex f2;
        if(edgeFaces.size() > 1)
            f2 = *(++edgeFaces.begin());

        // ignore edges inside of a cell
        if(mesh.getLabel(indexAttr[f1].label) == mesh.getLabel(indexAttr[f2].label))
            continue;
        //#pragma omp critical
        {
            wallAreas[std::make_pair(indexAttr[f1].label, indexAttr[f2].label)] += edgeLength;
            wallAreas[std::make_pair(indexAttr[f2].label, indexAttr[f1].label)] += edgeLength;
            wallEdges[std::make_pair(indexAttr[f1].label, indexAttr[f2].label)].insert(e);
            wallEdges[std::make_pair(indexAttr[f2].label, indexAttr[f1].label)].insert(e);
        }
    }
}


bool Tissue::initialize(bool create_dual, bool extended_dual) {

    if(parm("Verbose") == "True")
        mdxInfo << "Tissue::initialize" << endl;
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh or mesh->file().isEmpty())
        throw(QString("Tissue::initialize No current mesh"));

    // Get the tissue and cell graph names
    int cellDim = parm("Dimension").toInt();
    TissueName = parm("Tissue");
    TissueDualName = parm("Tissue Dual");

    cellTissue.processTissueParms(*this);
    if(mesh->ccStructure(TissueName).maxDimension != cellDim)
        mesh->ccStructure(TissueName) = CCStructure(cellDim);
    if(mesh->ccStructure(TissueDualName).maxDimension < 2)
        mesh->ccStructure(TissueDualName) = CCStructure(2);

    cellTissue.initialize(
        &mesh->ccStructure(TissueName), &mesh->ccStructure(TissueDualName), &mesh->indexAttr());

    // necessary for createdualextended
    restore(mesh->ccStructure(TissueName));

    // restore the tissue if requested
    if(create_dual) {
        if(parm("Verbose") == "True")
            mdxInfo << "Tissue::initialize: Dual graph restored" << endl;
        if(!extended_dual) {
            cellTissue.createDual();
            cellTissue.updateGeometry();
        }
        else
            createDualExtended( mesh->ccStructure(TissueName),  mesh->ccStructure(TissueDualName));
    }

    // Setup the complex attributes
    mesh->ccAttr(TissueName, "Type") = "Tissue";
    mesh->ccAttr(TissueName, "Dual") = TissueDualName;
    mesh->ccAttr(TissueDualName, "Type") = "TissueDual";
    mesh->ccAttr(TissueDualName, "Dual") = TissueName;

    mesh->updateAll();


    return true;
}

// Rewritten function from the original
void Tissue::createDualExtended(CCStructure &cs, CCStructure &csDual) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("Tissue::createDualExtended No current mesh"));

    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");
    csDual.clear();
    for(uint i = 0; i < cellAttr.size(); i++) {
        auto it = cellAttr.begin();
        advance(it, i);
        Tissue::CellData& cD = it->second;
        CCIndex v = CCIndexFactory.getIndex();
        csDual.addCell(v);
        indexAttr[v].pos = cD.centroid;
        vMAttr[v].dualCell = &cD;
        cD.dualVertex = v;
    }
    for(uint i = 0; i < cellAttr.size(); i++) {
        auto it = cellAttr.begin();
        advance(it, i);
        Tissue::CellData& cD = it->second;
        std::set<int> labels;
        for(CCIndex e : cD.perimeterEdges)
            for(CCIndex f : cs.incidentCells(e, 2))
                if(indexAttr[f].label != cD.label)
                    labels.insert(indexAttr[f].label);
        for(int label : labels) {
            CCIndex e = CCIndexFactory.getIndex();
            csDual.addCell(e, +cD.dualVertex -cellAttr[label].dualVertex);
        }
    }
}

// basically restore tissue's properties
void Tissue::restore(CCStructure csCurr) {

    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("CellTissueProcess::restore No current mesh"));
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
    Tissue::FaceDataAttr& faceAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");

    if(parm("Verbose") == "True")
        mdxInfo << "Tissue::restore" << endl;

    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    CCStructure& cs = mesh->ccStructure("Tissue");
    CCStructure& csDual = mesh->ccStructure("TissueDual");
    CCStructure& csVisual = mesh->ccStructure("TissueVisual");
    csVisual.clear();
    updateGeometry(cs, indexAttr);
    updateGeometry(csDual, indexAttr);
    updateGeometry(csVisual, indexAttr);

    // perimeters and borders lengths between cells
    neighborhood2D(*mesh, cs, indexAttr, wallAreas, wallEdges);

    // initialize cellAttr, at the beginning, owners are based on labels
    std::set<int> labels = getAllLabelsFaces(cs, indexAttr);
    for(int label : labels) {
        if(label < 1)
            throw(QString("Tissue::restore: cell with unallowed label " + label));
        // if a cell with the same label exist, we just reset its faces set,
        // which will be refilled in the next loop, otherwise, create a new cell
        if(cellAttr.find(label) != cellAttr.end()) {
            // verify that the cell was actually initialized
            if(!cellAttr[label].cellFaces)
                cellAttr[label].cellFaces = new std::set<CCIndex>();
        } else {
            cellAttr[label] = Tissue::CellData();
            cellAttr[label].tissue = this;
            cellAttr[label].label = label;
            if(parm("Verbose") == "True")
                mdxInfo << "Tissue::restore: A new cell has been created, label: " << label << endl;
        }
    }

    // cells that are not present anymore should be removed
    std::set<int> to_delete;
    for(auto c : cellAttr) {
        if(labels.find(c.second.label) == labels.end() || c.second.label < 1)
            to_delete.insert(c.second.label);
    }
    for(int label : to_delete) {
        Tissue::CellType type = cellAttr[label].type; /// IMPORTANT, if you refer to cellAttr the cell is recreated!
        cellAttr.erase(label);
        if(parm("Verbose") == "True")
            mdxInfo << "Tissue::restore: cell " << label << ", " << Tissue::ToString(type) <<  " not in the system anymore and therefore removed" << endl;
    }

    // restore CELLS properties
    for(uint i = 0; i < cellAttr.size(); i++) {
        auto it = cellAttr.begin();
        advance(it, i);
        Tissue::CellData& cD = it->second;
        cD.restore(mesh, this, cs, csDual, indexAttr, faceAttr, edgeAttr, vMAttr);
        cD.update(cs, indexAttr, faceAttr, edgeAttr, vMAttr, EPS);
        // restore the MF visuals
        if(parm("Draw MF") == "True") {
            cD.a1_v1 = CCIndexFactory.getIndex();
            csVisual.addCell(cD.a1_v1);
            cD.a1_v2 = CCIndexFactory.getIndex();
            csVisual.addCell(cD.a1_v2);
            cD.a1_e = CCIndexFactory.getIndex();
            csVisual.addCell(cD.a1_e, +cD.a1_v1 - cD.a1_v2);
            cD.a2_v1 = CCIndexFactory.getIndex();
            csVisual.addCell(cD.a2_v1);
            cD.a2_v2 = CCIndexFactory.getIndex();
            csVisual.addCell(cD.a2_v2);
            cD.a2_e = CCIndexFactory.getIndex();
            csVisual.addCell(cD.a2_e, +cD.a2_v1 - cD.a2_v2);
        }
       if(parm("Draw Bounding Box") == "True") {
            cD.box_v[0] = CCIndexFactory.getIndex();
            cD.box_v[1] = CCIndexFactory.getIndex();
            cD.box_v[2] = CCIndexFactory.getIndex();
            cD.box_v[3] = CCIndexFactory.getIndex();
            cD.box_e[0] = CCIndexFactory.getIndex();
            cD.box_e[1] = CCIndexFactory.getIndex();
            cD.box_e[2] = CCIndexFactory.getIndex();
            cD.box_e[3] = CCIndexFactory.getIndex();
            cD.axisCenter_v = CCIndexFactory.getIndex();
            cD.axisMax_v = CCIndexFactory.getIndex();
            cD.axisMin_v = CCIndexFactory.getIndex();
            cD.axisMax_e = CCIndexFactory.getIndex();
            cD.axisMin_e = CCIndexFactory.getIndex();
            csVisual.addCell(cD.box_v[0]);
            csVisual.addCell(cD.box_v[1]);
            csVisual.addCell(cD.box_v[2]);
            csVisual.addCell(cD.box_v[3]);
            csVisual.addCell(cD.box_e[0], +cD.box_v[0] - cD.box_v[1]);
            csVisual.addCell(cD.box_e[1], +cD.box_v[1] - cD.box_v[2]);
            csVisual.addCell(cD.box_e[2], +cD.box_v[2] - cD.box_v[3]);
            csVisual.addCell(cD.box_e[3], +cD.box_v[3] - cD.box_v[0]);
            csVisual.addCell(cD.axisCenter_v);
            csVisual.addCell(cD.axisMax_v);
            csVisual.addCell(cD.axisMin_v);
            csVisual.addCell(cD.axisMax_e, +cD.axisCenter_v - cD.axisMax_v);
            csVisual.addCell(cD.axisMin_e, +cD.axisCenter_v - cD.axisMin_v);
        }
        if(parm("Draw Auxin-Flux Rescale").toDouble() > 0) {
            cD.auxinFlux_v1 = CCIndexFactory.getIndex();
            csVisual.addCell(cD.auxinFlux_v1);
            cD.auxinFlux_v2 = CCIndexFactory.getIndex();
            csVisual.addCell(cD.auxinFlux_v2);
            cD.auxinFlux_e = CCIndexFactory.getIndex();
            csVisual.addCell(cD.auxinFlux_e, +cD.auxinFlux_v1 - cD.auxinFlux_v2);
            cD.auxinFlux_v3 = CCIndexFactory.getIndex();
            csVisual.addCell(cD.auxinFlux_v3);
            cD.auxinFlux_v4 = CCIndexFactory.getIndex();
            csVisual.addCell(cD.auxinFlux_v4);
            cD.auxinFlux_f = CCIndexFactory.getIndex();
            std::vector<CCIndex> vs;
            vs.push_back(cD.auxinFlux_v2);
            vs.push_back(cD.auxinFlux_v3);
            vs.push_back(cD.auxinFlux_v4);
            addFace(csVisual, cD.auxinFlux_f,  vs);
            indexAttr[cD.auxinFlux_f].label = 1;
        }
        if(parm("Draw PGD Rescale").toDouble() > 0) {
            cD.PDGmax_v1 = CCIndexFactory.getIndex();
            cD.PDGmax_v2 = CCIndexFactory.getIndex();
            cD.PDGmin_v1 = CCIndexFactory.getIndex();
            cD.PDGmin_v2 = CCIndexFactory.getIndex();
            csVisual.addCell(cD.PDGmax_v1);
            csVisual.addCell(cD.PDGmax_v2);
            csVisual.addCell(cD.PDGmin_v1);
            csVisual.addCell(cD.PDGmin_v2);
            cD.PDGmax_e = CCIndexFactory.getIndex();
            csVisual.addCell(cD.PDGmax_e, +cD.PDGmax_v1 - cD.PDGmax_v2);
            cD.PDGmin_e = CCIndexFactory.getIndex();
            csVisual.addCell(cD.PDGmin_e, +cD.PDGmin_v1 - cD.PDGmin_v2);
        }
    }

    // restore FACES previous properties and initialize them
    for(uint i = 0; i < cs.faces().size(); i++) {
        CCIndex f = cs.faces()[i];
        Tissue::FaceData& fD = faceAttr[f];
        fD.restore(f, cs, indexAttr);
        fD.update(f, cs, indexAttr, vMAttr);
        if(cellAttr[indexAttr[f].label].selected)
            indexAttr[f].selected = true;
        // restore the tensor strain rate visuals
        fD.G_v0 = CCIndexFactory.getIndex();
        csVisual.addCell(fD.G_v0);
        fD.G_v1 = CCIndexFactory.getIndex();
        csVisual.addCell(fD.G_v1);
        fD.G_v2 = CCIndexFactory.getIndex();
        csVisual.addCell(fD.G_v2);
        fD.G_e1 = CCIndexFactory.getIndex();
        csVisual.addCell(fD.G_e1, +fD.G_v0 - fD.G_v1);
        fD.G_e2 = CCIndexFactory.getIndex();
        csVisual.addCell(fD.G_e2, +fD.G_v0 - fD.G_v2);
        if(parm("Draw MF") == "True") {
            fD.a1_v1 = CCIndexFactory.getIndex();
            csVisual.addCell(fD.a1_v1);
            fD.a1_v2 = CCIndexFactory.getIndex();
            csVisual.addCell(fD.a1_v2);
            fD.a1_e = CCIndexFactory.getIndex();
            csVisual.addCell(fD.a1_e, +fD.a1_v1 - fD.a1_v2);
            fD.a2_v1 = CCIndexFactory.getIndex();
            csVisual.addCell(fD.a2_v1);
            fD.a2_v2 = CCIndexFactory.getIndex();
            csVisual.addCell(fD.a2_v2);
            fD.a2_e = CCIndexFactory.getIndex();
            csVisual.addCell(fD.a2_e);
            fD.a2_e = CCIndexFactory.getIndex();
            csVisual.addCell(fD.a2_e, +fD.a2_v1 - fD.a2_v2);

        }
    }

    // Restore EDGES
    for(uint i = 0; i < cs.edges().size(); i++) {
        CCIndex e = cs.edges()[i];
        Tissue::EdgeData& eD = edgeAttr[e];
        eD.restore(e, cs, indexAttr, cellAttr);
        eD.update(e, cs, indexAttr, vMAttr, faceAttr);
        if(parm("Draw PINs").toDouble() > 0) {
            if(eD.type == Wall) {
                    eD.pin_v1_1 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_v1_1);
                    eD.pin_v2_1 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_v2_1);
                    eD.pin_v3_1 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_v3_1);
                    eD.pin_v4_1 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_v4_1);
                    eD.pin_e1_1 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_e1_1, + eD.pin_v1_1 - eD.pin_v2_1);
                    eD.pin_e2_1 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_e2_1, + eD.pin_v2_1 - eD.pin_v3_1);
                    eD.pin_e3_1 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_e3_1, + eD.pin_v3_1 - eD.pin_v4_1);
                    eD.pin_e4_1 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_e4_1, + eD.pin_v4_1 - eD.pin_v1_1);
                    std::vector<CCIndex> vs;
                    vs.push_back(eD.pin_v1_1);
                    vs.push_back(eD.pin_v2_1);
                    vs.push_back(eD.pin_v3_1);
                    vs.push_back(eD.pin_v4_1);
                    eD.pin_f_1 = CCIndexFactory.getIndex();
                    addFace(csVisual, eD.pin_f_1,  vs);
                    indexAttr[eD.pin_f_1].label = 2; // to color pins

                    eD.pin_v1_2 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_v1_2);
                    eD.pin_v2_2 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_v2_2);
                    eD.pin_v3_2 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_v3_2);
                    eD.pin_v4_2 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_v4_2);
                    eD.pin_e1_2 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_e1_2, + eD.pin_v1_2 - eD.pin_v2_2);
                    eD.pin_e2_2 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_e2_2, + eD.pin_v2_2 - eD.pin_v3_2);
                    eD.pin_e3_2 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_e3_2, + eD.pin_v3_2 - eD.pin_v4_2);
                    eD.pin_e4_2 = CCIndexFactory.getIndex();
                    csVisual.addCell(eD.pin_e4_2, + eD.pin_v4_2 - eD.pin_v1_2);
                    vs.clear();
                    vs.push_back(eD.pin_v1_2);
                    vs.push_back(eD.pin_v2_2);
                    vs.push_back(eD.pin_v3_2);
                    vs.push_back(eD.pin_v4_2);
                    eD.pin_f_2 = CCIndexFactory.getIndex();
                    addFace(csVisual, eD.pin_f_2,  vs);
                    indexAttr[eD.pin_f_2].label = 2; // to color pins


            }
        }
    }


    // restore VERTICES
    for(uint i = 0; i < cs.vertices().size(); i++) {
        CCIndex v = cs.vertices()[i];
        Tissue::VertexData& vD = vMAttr[v];
        vD.restore(v, cs, indexAttr, cellAttr, edgeAttr);
        vMAttr[v].prevPos = indexAttr[v].pos;
        vD.V_v0 = CCIndexFactory.getIndex();
        csVisual.addCell(vD.V_v0);
        vD.V_v1 = CCIndexFactory.getIndex();
        csVisual.addCell(vD.V_v1);
        vD.V_e = CCIndexFactory.getIndex();
        csVisual.addCell(vD.V_e, +vD.V_v0 - vD.V_v1);
    }

    // clear unused edges
    std::set<CCIndex> old_cells;
    for(auto i : edgeAttr)
        if(!cs.hasCell(i.first) && !csVisual.hasCell(i.first) && !csDual.hasCell(i.first))
            old_cells.insert(i.first);
    for(CCIndex e : old_cells)
            edgeAttr.erase(e);
    old_cells.clear();
    // clear unused vertices
    for(auto i : vMAttr)
        if(!cs.hasCell(i.first) && !csVisual.hasCell(i.first) && !csDual.hasCell(i.first))
            old_cells.insert(i.first);
    for(CCIndex v : old_cells)
            vMAttr.erase(v);
    old_cells.clear();
    // clear unused faces
    for(auto i : faceAttr)
        if(!cs.hasCell(i.first) && !csVisual.hasCell(i.first) && !csDual.hasCell(i.first))
            old_cells.insert(i.first);
    for(CCIndex f : old_cells)
            faceAttr.erase(f);
    old_cells.clear();
    // clear unused indexes
    for(auto i : indexAttr)
        if(!cs.hasCell(i.first) && !csVisual.hasCell(i.first) && !csDual.hasCell(i.first)) // IMPORTANT: you have to check all the cell complexes
            old_cells.insert(i.first);
    for(CCIndex i : old_cells)
            indexAttr.erase(i);

    mesh->updateAll();
}

// update geometries, signals....
bool Tissue::step(double Dt) {

    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("CellTissueProcess::initialize No current mesh"));

    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
    Tissue::FaceDataAttr& faceAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");

    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    CCStructure& cs = mesh->ccStructure("Tissue");
    CCStructure& csDual = mesh->ccStructure("TissueDual");
    CCStructure& csVisual = mesh->ccStructure("TissueVisual");
    updateGeometry(cs, indexAttr);
    updateGeometry(csDual, indexAttr);
    updateGeometry(csVisual, indexAttr);


    CCIndexIntAttr& cellTypeSignal = mesh->signalAttr<int>("Cell Type");
    CCIndexIntAttr& cellLineageSignal = mesh->signalAttr<int>("Cell Lineage");
    CCIndexDoubleAttr& divisionCountSignal = mesh->signalAttr<double>("Division Count");
    CCIndexDoubleAttr& auxinSignal = mesh->signalAttr<double>("Chems: Auxin");
    CCIndexDoubleAttr& auxinByAreaSignal = mesh->signalAttr<double>("Chems: Auxin By Area");
    CCIndexDoubleAttr& intercellularAuxinSignal = mesh->signalAttr<double>("Chems: Intercellular Auxin");
    CCIndexDoubleAttr& Aux1CytSignal = mesh->signalAttr<double>("Chems: Aux1 Cyt");
    CCIndexDoubleAttr& Aux1MemSignal = mesh->signalAttr<double>("Chems: Aux1 Mem");
    CCIndexDoubleAttr& MFImpactSignal = mesh->signalAttr<double>("Chems: MF Impact");
    CCIndexDoubleAttr& geomImpactSignal = mesh->signalAttr<double>("Chems: Geometry Impact");
    CCIndexDoubleAttr& auxinRatioSignal = mesh->signalAttr<double>("Chems: Auxin Transport");
    CCIndexDoubleAttr& auxinGradSignal = mesh->signalAttr<double>("Chems: Auxin Gradient");
    CCIndexDoubleAttr& auxinFluxImpactSignal = mesh->signalAttr<double>("Chems: Auxin Flux Impact");
    CCIndexDoubleAttr& divInhibitorCytSignal = mesh->signalAttr<double>("Chems: Division Inhibitor by Area");
    CCIndexDoubleAttr& divPromoterCytSignal = mesh->signalAttr<double>("Chems: Division Promoter by Area");
    CCIndexDoubleAttr& divProbSignal = mesh->signalAttr<double>("Chems: Division Probability");
    CCIndexDoubleAttr& divTimeSignal = mesh->signalAttr<double>("Chems: Division Time");
    CCIndexDoubleAttr& Pin1CytSignal = mesh->signalAttr<double>("Chems: Pin1 Cyt");
    CCIndexDoubleAttr& Pin1MemSignal = mesh->signalAttr<double>("Chems: Pin1 Mem");
    CCIndexDoubleAttr& pin1SensitivitySignal = mesh->signalAttr<double>("Chems: Pin Sensitivity");
    CCIndexDoubleAttr& PINOIDCytSignal = mesh->signalAttr<double>("Chems: PINOID Cyt");
    CCIndexDoubleAttr& PP2ACytSignal = mesh->signalAttr<double>("Chems: PP2A Cyt");
    CCIndexDoubleAttr& PINOIDMemSignal = mesh->signalAttr<double>("Chems: PINOID Mem");
    CCIndexDoubleAttr& PP2AMemSignal = mesh->signalAttr<double>("Chems: PP2A Mem");
    CCIndexDoubleAttr& pressureSignal = mesh->signalAttr<double>("Mechs: Turgor Pressure");
    CCIndexDoubleAttr& sigmaASignal = mesh->signalAttr<double>("Mechs: Anisotropic Force");
    CCIndexDoubleAttr& edgeStiffnessSignal = mesh->signalAttr<double>("Mechs: Edge Stiffness");
    CCIndexDoubleAttr& edgeStrainSignal = mesh->signalAttr<double>("Mechs: Edge Strain Rate");
    CCIndexDoubleAttr& growthRateSignal = mesh->signalAttr<double>("Mechs: Growth Rate");
    CCIndexDoubleAttr& distanceConstrainSignal = mesh->signalAttr<double>("PBD: distance constrain");
    CCIndexDoubleAttr& bendingConstrainSignal = mesh->signalAttr<double>("PBD: bending constrain");
    CCIndexDoubleAttr& shapeConstrainSignal = mesh->signalAttr<double>("PBD: shape constrain");
    CCIndexDoubleAttr& strainConstrainSignal = mesh->signalAttr<double>("PBD: strain constrain");
    CCIndexDoubleAttr& pressureConstrainSignal = mesh->signalAttr<double>("PBD: pressure constrain");

    cellTypeSignal.clear();
    cellLineageSignal.clear();
    divisionCountSignal.clear();
    auxinSignal.clear();
    auxinByAreaSignal.clear();
    intercellularAuxinSignal.clear();
    Aux1CytSignal.clear();
    Aux1MemSignal.clear();
    MFImpactSignal.clear();
    geomImpactSignal.clear();
    auxinFluxImpactSignal.clear();
    auxinRatioSignal.clear();
    auxinGradSignal.clear();
    divInhibitorCytSignal.clear();
    divPromoterCytSignal.clear();
    divProbSignal.clear();
    divTimeSignal.clear();
    Pin1CytSignal.clear();
    Pin1MemSignal.clear();
    pin1SensitivitySignal.clear();
    PINOIDCytSignal.clear();
    PP2ACytSignal.clear();
    PINOIDMemSignal.clear();
    PP2AMemSignal.clear();
    pressureSignal.clear();
    sigmaASignal.clear();
    edgeStrainSignal.clear();
    edgeStiffnessSignal.clear();
    growthRateSignal.clear();

    //updateGeometry(cs, indexAttr);
    //updateGeometry(csDual, indexAttr);


    // update TissueDual positions
    for(CCIndex v : csDual.vertices()) {
        Tissue::VertexData& vD = vMAttr[v];
        if(vD.dualCell)
            indexAttr[v].pos = ((CellData*)vD.dualCell)->centroid;
    }


    // update EDGES geometry
    //#pragma omp parallel for
    for(uint i = 0; i < cs.edges().size(); i++) {
        CCIndex e = cs.edges()[i];
        Tissue::EdgeData& eD = edgeAttr[e];
        eD.update(e, cs, indexAttr, vMAttr, faceAttr);
        // PINs
        double shift = 0.1; // distance from the walls
        if(parm("Draw PINs").toDouble() > 0) {
            double extension = parm("Draw PINs").toDouble();
            if(eD.type == Wall) {
                std::vector<CCIndex> fs;
                for(auto f : cs.incidentCells(e, 2))
                    fs.push_back(f);

                Tissue::FaceData fD = faceAttr[fs[0]];
                Point3d c = eD.outwardNormal[fs[0]];
                int label = fD.owner->label;
                double pin = eD.Pin1[label] + EPS;
                auto eb = cs.edgeBounds(e);
                CCIndex v1 = eb.first;
                CCIndex v2 = eb.second;
                Point3d versor = indexAttr[v2].pos - indexAttr[v1].pos;
                versor /= norm(versor);
                Point3d p1 = indexAttr[v1].pos + -c * shift + versor * shift;
                Point3d p2 = indexAttr[v2].pos + -c * shift - versor * shift;
                Point3d p4 = p1 + (-c * extension * pin) + versor * shift*3;
                Point3d p3 = p2 + (-c * extension  * pin) - versor * shift*3;
                indexAttr[eD.pin_v1_1].pos = p1+ Point3d(0,0,0.04);
                indexAttr[eD.pin_v2_1].pos = p2+ Point3d(0,0,0.04);
                indexAttr[eD.pin_v3_1].pos = p3+ Point3d(0,0,0.04);
                indexAttr[eD.pin_v4_1].pos = p4+ Point3d(0,0,0.04);

                if(fs.size() > 1) {
                    fD = faceAttr[fs[1]];
                    c = eD.outwardNormal[fs[1]];
                    label = fD.owner->label;
                    pin = eD.Pin1[label] + EPS;
                    p1 = indexAttr[v1].pos + -c * shift + versor * shift;
                    p2 = indexAttr[v2].pos + -c * shift - versor * shift;
                    p4 = p1 + (-c * extension * pin) + versor * shift*3;
                    p3 = p2 + (-c * extension * pin) - versor * shift*3;
                    indexAttr[eD.pin_v1_2].pos = p1 + Point3d(0,0,0.04);
                    indexAttr[eD.pin_v2_2].pos = p2 + Point3d(0,0,0.04);
                    indexAttr[eD.pin_v3_2].pos = p3 + Point3d(0,0,0.04);
                    indexAttr[eD.pin_v4_2].pos = p4 + Point3d(0,0,0.04);
                }
            }
        }
     }

    // update VERTICES, forces etc...visuals
    //#pragma omp parallel for
    for(uint i = 0; i < cs.vertices().size(); i++) {
        CCIndex v = cs.vertices()[i];
        Tissue::VertexData& vD = vMAttr[v];
        Point3d totalForce, sigmaForce;
        distanceConstrainSignal[v] = norm(vD.corrections["distance"]);
        bendingConstrainSignal[v] = norm(vD.corrections["bending"]);
        shapeConstrainSignal[v] = norm(vD.corrections["shape"]);
        strainConstrainSignal[v] = norm(vD.corrections["strain"]);
        pressureConstrainSignal[v] = norm(vD.corrections["pressure"]);
        for(auto m : vD.forces) {
            if(std::get<1>(m) == QString("sigmaAY"))
                sigmaForce += std::get<2>(m);
            totalForce += std::get<2>(m);
        }
        if(parm("Draw Velocity Rescale").toDouble() > 0) {
            indexAttr[vD.V_v0].pos = indexAttr[v].pos;
            indexAttr[vD.V_v1].pos = indexAttr[vD.V_v0].pos + vD.velocity * parm("Draw Velocity Rescale").toDouble();
        }
    }

    // update basic CELL data and visuals (has to go before faces as we are setting cD.selected, this it limit iterations)
    for(uint i = 0; i < cellAttr.size(); i++) {
        auto it = cellAttr.begin();
        advance(it, i);
        Tissue::CellData& cD = it->second;
        cD.update(cs, indexAttr, faceAttr, edgeAttr, vMAttr, Dt);
        cD.selected = false;
        if(parm("Draw MF") == "True") {
            indexAttr[cD.a1_v1].pos = cD.centroid;
            indexAttr[cD.a1_v2].pos = cD.centroid + cD.a1;
            indexAttr[cD.a2_v1].pos = cD.centroid;
            indexAttr[cD.a2_v2].pos = cD.centroid + cD.a2;
        }
        if(parm("Draw Auxin-Flux Rescale").toDouble() > 0) { // Draw auxin flow arrows
            indexAttr[cD.auxinFlux_v1].pos = cD.centroid;
            indexAttr[cD.auxinFlux_v2].pos = cD.centroid + cD.auxinFluxVector  * parm("Draw Auxin-Flux Rescale").toDouble();
            indexAttr[cD.auxinFlux_v3].pos = indexAttr[cD.auxinFlux_v2].pos;
            indexAttr[cD.auxinFlux_v3].pos = mdx::rotatePoint2D(indexAttr[cD.auxinFlux_v3].pos, cD.centroid, 0.15);
            Point3d v = cD.centroid - indexAttr[cD.auxinFlux_v3].pos;
            indexAttr[cD.auxinFlux_v3].pos += v / norm(v) * (0.35 * norm(v));;
            indexAttr[cD.auxinFlux_v4].pos = indexAttr[cD.auxinFlux_v2].pos;
            indexAttr[cD.auxinFlux_v4].pos = mdx::rotatePoint2D(indexAttr[cD.auxinFlux_v4].pos, cD.centroid, -0.15);
            v = cD.centroid - indexAttr[cD.auxinFlux_v4].pos;
            indexAttr[cD.auxinFlux_v4].pos += v / norm(v) * (0.35 * norm(v));
            indexAttr[cD.auxinFlux_v1].pos[2] += 0.1;
            indexAttr[cD.auxinFlux_v2].pos[2] += 0.1;
            indexAttr[cD.auxinFlux_v3].pos[2] += 0.2;
            indexAttr[cD.auxinFlux_v4].pos[2] += 0.2;
        }
        if(parm("Draw Bounding Box") == "True") {
            indexAttr[cD.box_v[0]].pos = cD.box[0];
            indexAttr[cD.box_v[1]].pos = cD.box[1];
            indexAttr[cD.box_v[2]].pos = cD.box[2];
            indexAttr[cD.box_v[3]].pos = cD.box[3];
            indexAttr[cD.axisCenter_v].pos = cD.centroid;
            indexAttr[cD.axisMax_v].pos = cD.centroid + cD.axisMax;
            indexAttr[cD.axisMin_v].pos = cD.centroid + cD.axisMin;
        }
        if(parm("Draw PGD Rescale").toDouble() > 0) {
            double axixMax_gr = cD.axisMax_growthRate * parm("Draw PGD Rescale").toDouble();
            double axixMin_gr = cD.axisMin_growthRate * parm("Draw PGD Rescale").toDouble();
            indexAttr[cD.PDGmax_v1].pos = cD.centroid - (cD.axisMax / norm(cD.axisMax)) * axixMax_gr ;
            indexAttr[cD.PDGmax_v2].pos = cD.centroid + (cD.axisMax / norm(cD.axisMax)) * axixMax_gr ;
            indexAttr[cD.PDGmin_v1].pos = cD.centroid - (cD.axisMin / norm(cD.axisMin)) * axixMin_gr ;
            indexAttr[cD.PDGmin_v2].pos = cD.centroid + (cD.axisMin / norm(cD.axisMin)) * axixMin_gr ;
        }
    }

    // update FACES, MF, tensors etc...visuals
    //#pragma omp parallel for
    for(uint i = 0; i < cs.faces().size(); i++) {
        CCIndex f = cs.faces()[i];
        Tissue::FaceData& fD = faceAttr[f];
        Tissue::CellData& cD = cellAttr[indexAttr[f].label];
        fD.update(f, cs, indexAttr, vMAttr);
        if(indexAttr[f].selected)
            cD.selected = true;
        // Green Strain Tensor?
        if(parm("Draw Strain Tensors Rescale").toDouble() > 0) {
            indexAttr[fD.G_v0].pos = indexAttr[f].pos;
            indexAttr[fD.G_v1].pos = indexAttr[f].pos;
            indexAttr[fD.G_v2].pos = indexAttr[f].pos;
            indexAttr[fD.G_v1].pos += (Rotate(Point3d(cD.gMax,0,0), cD.gAngle) * parm("Draw Strain Tensors Rescale").toDouble());
            indexAttr[fD.G_v2].pos += (Rotate(Point3d(cD.gMin,0,0), cD.gAngle + M_PI/2) * parm("Draw Strain Tensors Rescale").toDouble());
        }
        // MF
        if(parm("Draw MF") == "True") {
            indexAttr[fD.a1_v1].pos = indexAttr[f].pos;
            indexAttr[fD.a1_v2].pos = indexAttr[f].pos + fD.a1;
            indexAttr[fD.a2_v1].pos = indexAttr[f].pos;
            indexAttr[fD.a2_v2].pos = indexAttr[f].pos + fD.a2;
        }

        // Signals
        cellTypeSignal[f] = cellAttr[indexAttr[f].label].type;
        cellLineageSignal[f] = cellAttr[indexAttr[f].label].lineage;
        divisionCountSignal[f] = cellAttr[indexAttr[f].label].divisionCount;
        auxinSignal[f] = cD.auxin;
        auxinByAreaSignal[f] = cD.auxin / cD.area;
        Aux1CytSignal[f] = fD.Aux1Cyt;
        Aux1MemSignal[f] = fD.Aux1Mem;
        divInhibitorCytSignal[f] = cD.divInhibitor / cD.area;
        divPromoterCytSignal[f] = cD.divPromoter / cD.area;
        divProbSignal[f] = cD.divProb;
        divTimeSignal[f] = exp(-0.1*cD.lastDivision);
        Pin1CytSignal[f] = cD.Pin1;
        Pin1MemSignal[f] = fD.Pin1Mem;
        PINOIDCytSignal[f] = cD.PINOID;
        PP2ACytSignal[f] = cD.PP2A;
        PINOIDMemSignal[f] = fD.PINOIDMem;
        PP2AMemSignal[f] = fD.PP2AMem;
        pressureSignal[f] = cD.pressure;
        //sigmaASignal[f] = norm(fD.sigmaA);
        sigmaASignal[f] = norm(fD.a1) - norm(fD.a2);

        fD.auxin = cD.auxin;
        fD.intercellularAuxin = 0;
        fD.Aux1Cyt = cD.Aux1;
        fD.Pin1Cyt = cD.Pin1;
        fD.auxinRatio = 0;
        fD.auxinGrad = 0;
        fD.pin1Sensitivity = 0;
        fD.Pin1Mem = 0;
        fD.Aux1Mem = 0;
        fD.PINOIDMem = 0;
        fD.PP2AMem = 0;
        fD.MFImpact = 0;
        fD.geomImpact = 0;
        fD.auxinFluxImpact = 0;
        fD.edgeStrain = 0;
        fD.edgeStiffness = 0;
        if(fD.type == Membrane) {
            int es = 0;
            for(CCIndex e : cs.incidentCells(f, 1)) {
                Tissue::EdgeData& eD = edgeAttr[e];
                if(eD.type == Tissue::Wall) {
                    fD.Pin1Mem += eD.Pin1[cD.label] / eD.length;
                    fD.Aux1Mem += eD.Aux1[cD.label] / eD.length;
                    fD.PINOIDMem += eD.PINOID[cD.label] / eD.length;
                    fD.PP2AMem += eD.PP2A[cD.label] / eD.length;
                    fD.intercellularAuxin += eD.intercellularAuxin / eD.length;
                    fD.edgeStrain += eD.strainRate;
                    fD.edgeStiffness += eD.eStiffness;
                    fD.pin1Sensitivity += eD.pin1Sensitivity[indexAttr[f].label];
                    fD.MFImpact += eD.MFImpact[indexAttr[f].label];
                    fD.geomImpact += eD.geomImpact[indexAttr[f].label];
                    fD.auxinFluxImpact += eD.auxinFluxImpact[indexAttr[f].label];
                    fD.auxinRatio += eD.auxinRatio[indexAttr[f].label];
                    fD.auxinGrad += eD.auxinGrad[indexAttr[f].label];
                    es++;
                }
            }
           fD.edgeStrain /= es;
           fD.edgeStiffness /= es;
           fD.pin1Sensitivity /= es;
           fD.MFImpact /= es;
           fD.geomImpact /= es;
           fD.auxinFluxImpact /= es;
           fD.auxinRatio /= es;
           fD.auxinGrad /= es;
           fD.Pin1Mem /= es;
           fD.Aux1Mem /= es;
           fD.PINOIDMem /= es;
           fD.PP2AMem /= es;
           fD.intercellularAuxin /= es;
        }
        intercellularAuxinSignal[f] = fD.intercellularAuxin;
        edgeStrainSignal[f] = fD.edgeStrain;
        edgeStiffnessSignal[f] = fD.edgeStiffness;
        pin1SensitivitySignal[f] = fD.pin1Sensitivity;
        MFImpactSignal[f] = fD.MFImpact;
        geomImpactSignal[f] = fD.geomImpact;
        auxinFluxImpactSignal[f] = fD.auxinFluxImpact;
        auxinGradSignal[f] = fD.auxinGrad;
        auxinRatioSignal[f] = fD.auxinRatio;
        growthRateSignal[f] = fD.growthRate;
    }

    // perimeters and borders lengths between cells
    neighborhood2D(*mesh, cs, indexAttr, wallAreas, wallEdges);

    return false;
}

void Tissue::CellData::division(const CCStructure &cs,
                                CellDataAttr& cellAttr, FaceDataAttr& faceAttr, EdgeDataAttr& edgeAttr,
                                CellData& cD1, CellData& cD2, std::map<CellType, int> maxAreas, bool ignoreCellType) {

    if(!tissue)
        throw(QString("Tissue::CellData::division: tissue not set"));

    Point3d cm, QCcm;
    for(auto c : cellAttr) {
        Tissue::CellData& cDn = cellAttr[c.first];
        cm += cDn.centroid;
        if(cDn.type == QC)
            QCcm += cDn.centroid;
    }
    cm /= cellAttr.size();
    QCcm /= 2; // assumes two QC cells
    Point3d rootAxisY0 = Point3d(QCcm[0], -BIG_VAL, 0);
    Point3d rootAxisY1 = Point3d(QCcm[0], BIG_VAL, 0);
    Point3d rootAxisX0 = Point3d(-BIG_VAL, QCcm[1], 0);
    Point3d rootAxisX1 = Point3d(BIG_VAL, QCcm[1], 0);
    CellType type1 = type, type2 = type;
    bool periclinalDivision1 = false, periclinalDivision2 = false;
    if(!ignoreCellType) {
        if(type == CEI) {
            double dist1 = DistancePtLine(rootAxisX0, rootAxisX1, cD1.centroid);
            double dist2 = DistancePtLine(rootAxisX0, rootAxisX1, cD2.centroid);
            if(dist1 > dist2) {
                type1 = CEID;
                type2 = CEI;
                periclinalDivision1 = true;
                periclinalDivision2 = false;
            } else {
                type2 = CEID;
                type1 = CEI;
                periclinalDivision1 = false;
                periclinalDivision2 = true;
            }
        } else if(type == VascularInitial) {
            double dist1 = DistancePtLine(rootAxisX0, rootAxisX1, cD1.centroid);
            double dist2 = DistancePtLine(rootAxisX0, rootAxisX1, cD2.centroid);
            if(dist1 > dist2) {
                type1 = Vascular;
                type2 = VascularInitial;
            } else {
                type2 = Vascular;
                type1 = VascularInitial;
            }
        } else if(type == CEID) {
            double dist1 = DistancePtLine(rootAxisY0, rootAxisY1, cD1.centroid);
            double dist2 = DistancePtLine(rootAxisY0, rootAxisY1, cD2.centroid);
            if(dist1 > dist2) {
                type1 = Cortex;
                type2 = Endodermis;
            } else {
                type2 = Cortex;
                type1 = Endodermis;
            }
        } else if(type == EpLrcInitial) {
            CellType nextType;
            if(periclinalDivision) {
                nextType = LRC;
                double dist1 = cD1.centroid.y();//norm(cD1.centroid - QCcm);
                double dist2 = cD2.centroid.y();//norm(cD2.centroid - QCcm);
                if(dist1 < dist2) {
                    type1 = nextType;
                    type2 = EpLrcInitial;
                    periclinalDivision1 = false;
                    periclinalDivision2 = !periclinalDivision;
                }else  {
                    type2 = nextType;
                    type1 = EpLrcInitial;
                    periclinalDivision1 = !periclinalDivision;
                    periclinalDivision2 = false;
                }
            }
            else {
                nextType = Epidermis;
                double dist1 = cD1.centroid.y();//norm(cD1.centroid - QCcm);
                double dist2 = cD2.centroid.y();//norm(cD2.centroid - QCcm);
                if(dist1 > dist2) {
                    type1 = nextType;
                    type2 = EpLrcInitial;
                    periclinalDivision1 = false;
                    periclinalDivision2 = !periclinalDivision;
                }else  {
                    type2 = nextType;
                    type1 = EpLrcInitial;
                    periclinalDivision1 = !periclinalDivision;
                    periclinalDivision2 = false;
                }
            }

        } else if(type == ColumellaInitial) {
            double dist1 = DistancePtLine(rootAxisX0, rootAxisX1, cD1.centroid);
            double dist2 = DistancePtLine(rootAxisX0, rootAxisX1, cD2.centroid);
            if(dist1 > dist2) {
                type1 = Columella;
                type2 = ColumellaInitial;
            } else {
                type2 = Columella;
                type1 = ColumellaInitial;
            }
        }
    }

    cD1.tissue = tissue;
    cD1.type = type1;
    cD1.lineage = lineage;
    cD1.periclinalDivision = periclinalDivision1;
    if(mfDelete && mfRORate > 0) {
        cD1.a1 = Point3d(0, 0.0001, 0);
        cD1.a2 = Point3d(0.0001, 0, 0);
    } else {
        cD1.a1 = a1;
        cD1.a2 = a2;
    }
    cD1.cellMaxArea = maxAreas[cD1.type];
    cD1.restArea = (cD1.area / area) * restArea;
    cD1.prevArea = (cD1.area / area) * prevArea;
    cD1.mfRORate = mfRORate;
    cD1.divAlg = divAlg;
    cD1.auxin = auxin / 2;
    cD1.auxinProdRate = auxinProdRate;
    cD1.pinProdRate = pinProdRate;
    cD1.pinInducedRate = pinInducedRate;
    cD1.aux1ProdRate = aux1ProdRate;
    cD1.aux1InducedRate = aux1InducedRate;
    cD1.aux1MaxEdge = aux1MaxEdge;
    cD1.growthFactor = growthFactor;
    cD1.pressureMax = pressureMax;
    cD1.Pin1 = Pin1 / 2;
    cD1.Aux1 = Aux1 / 2;
    cD1.PINOID = PINOID / 2;
    cD1.PP2A = PP2A / 2;
    cD1.divPromoter = divPromoter / 2;
    cD1.divisionCount = divisionCount+1;

    cD2.tissue = tissue;
    cD2.type = type2;
    cD2.lineage = lineage;
    cD2.periclinalDivision = periclinalDivision2;
    if(mfDelete && mfRORate > 0) {
        cD2.a1 = Point3d(0, 0.0001, 0);
        cD2.a2 = Point3d(0.0001, 0, 0);
    } else {
        cD2.a1 = a1;
        cD2.a2 = a2;
    }
    cD2.cellMaxArea = maxAreas[cD2.type];
    cD2.restArea = (cD2.area / area) * restArea;
    cD2.prevArea = (cD2.area / area) * prevArea;
    cD2.mfRORate = mfRORate;
    cD2.divAlg = divAlg;
    cD2.auxin = auxin / 2;
    cD2.auxinProdRate = auxinProdRate;
    cD2.pinProdRate = pinProdRate;
    cD2.pinInducedRate = pinInducedRate;
    cD2.aux1ProdRate = aux1ProdRate;
    cD2.aux1InducedRate = aux1InducedRate;
    cD2.aux1MaxEdge = aux1MaxEdge;
    cD2.growthFactor = growthFactor;
    cD2.pressureMax = pressureMax;
    cD2.Pin1 = Pin1 / 2;
    cD2.Aux1 = Aux1 / 2;
    cD2.PINOID = PINOID / 2;
    cD2.PP2A = PP2A / 2;
    cD2.divInhibitor = 0; //divInhibitor / 2;
    cD2.divPromoter = divPromoter / 2;
    cD2.divisionCount = divisionCount+1;

    std::vector<CCIndex> edges;
    edges.insert(edges.end(), cD1.perimeterEdges.begin(), cD1.perimeterEdges.end());
    edges.insert(edges.end(), cD2.perimeterEdges.begin(), cD2.perimeterEdges.end());
    for(CCIndex e : edges) {
        Tissue::EdgeData& eD = edgeAttr[e];
        eD.Pin1.erase(label);
        eD.auxinRatio.erase(label);
        eD.auxinGrad.erase(label);
        eD.Aux1.erase(label);
        eD.PINOID.erase(label);
        eD.PP2A.erase(label);
	eD.auxinFluxImpact.erase(label);
        eD.MFImpact.erase(label);
        eD.geomImpact.erase(label);
        eD.pin1SensitivityRaw.erase(label);
        eD.pin1Sensitivity.erase(label);
    }
    /*
    if(parm("Verbose") == "True")
        mdxInfo << "I, label " << label << " am a " << Tissue::ToString(type)
         << " and my daughters are " << Tissue::ToString(cD1.type) << " and " << Tissue::ToString(cD2.type)
         << " of max cell size " << cD1.cellMaxArea << " and " << cD2.cellMaxArea <<  endl;*/
}

