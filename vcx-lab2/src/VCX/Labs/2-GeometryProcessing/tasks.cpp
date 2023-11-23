#include <iostream>
#include <list>
#include <map>
#include <set>
#include <unordered_set>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtx/vector_angle.hpp>

#include <spdlog/spdlog.h>

#include "Labs/2-GeometryProcessing/DCEL.hpp"
#include "Labs/2-GeometryProcessing/tasks.h"

namespace VCX::Labs::GeometryProcessing {

#include "Labs/2-GeometryProcessing/marching_cubes_table.h"

    /******************* 1. Mesh Subdivision *****************/
    void SubdivisionMesh(Engine::SurfaceMesh const & input, Engine::SurfaceMesh & output, std::uint32_t numIterations) {
        // your code here
        output = input;
        for (uint32_t n = 0; n < numIterations; ++n) {
            DCEL links;
            links.AddFaces(output.Indices);
            if (! links.IsValid()) {
                printf("1 : invalid");
                return;
            }
            std::vector<glm::vec3> pre_pos = output.Positions;
            output.Indices.clear();
            uint32_t                                                               length = output.Positions.size();
            std::map<std::pair<DCEL::VertexIdx, DCEL::VertexIdx>, DCEL::VertexIdx> m;
            for (DCEL::Triangle const & f : links.GetFaces()) {
                DCEL::VertexIdx index[6];
                for (int i = 0; i < 3; ++i) {
                    index[i]                               = f.Edges(i)->To();
                    DCEL::Vertex                 v         = links.GetVertex(index[i]);
                    std::vector<DCEL::VertexIdx> neighbors = v.GetNeighbors();
                    size_t                       len       = neighbors.size();
                    float                        u         = 0;
                    if (len == 3) u = 3.0 / 16;
                    else u = 3.0 / (8 * len);
                    glm::vec3 pos = pre_pos[index[i]];
                    pos           = pos * (1 - len * u);
                    for (size_t j = 0; j < len; ++j)
                        pos = pos + u * pre_pos[neighbors[j]];
                    output.Positions[index[i]] = pos;
                }
                for (int i = 0; i < 3; ++i) {
                    DCEL::HalfEdge const * e  = f.Edges(i);
                    DCEL::VertexIdx        v1 = e->OppositeVertex(), v2 = e->PairEdge()->OppositeVertex(), v3 = e->From(), v4 = e->To();
                    if (v3 > v4) {
                        DCEL::VertexIdx tmp = v3;
                        v3                  = v4;
                        v4                  = tmp;
                    }
                    auto p = std::make_pair(v3, v4);
                    if (m.find(p) != m.end()) {
                        index[3 + i] = m[p];
                    } else {
                        glm::vec3 pos = (float) 1 / 8 * (pre_pos[v1] + pre_pos[v2]) + (float) 3 / 8 * (pre_pos[v3] + pre_pos[v4]);
                        output.Positions.push_back(pos);
                        index[3 + i] = length;
                        m[p]         = length++;
                    }
                }
                output.Indices.push_back(index[0]);
                output.Indices.push_back(index[4]);
                output.Indices.push_back(index[3]);
                output.Indices.push_back(index[4]);
                output.Indices.push_back(index[1]);
                output.Indices.push_back(index[5]);
                output.Indices.push_back(index[3]);
                output.Indices.push_back(index[5]);
                output.Indices.push_back(index[2]);
                output.Indices.push_back(index[3]);
                output.Indices.push_back(index[4]);
                output.Indices.push_back(index[5]);
            }
        }
        return;
    }

    /******************* 2. Mesh Parameterization *****************/
    void Parameterization(Engine::SurfaceMesh const & input, Engine::SurfaceMesh & output, const std::uint32_t numIterations) {
        // your code here
        output = input;
        DCEL links;
        links.AddFaces(output.Indices);
        if (! links.IsValid()) {
            printf("2 : invalid");
            return;
        }
        // boundary and inside
        std::vector<DCEL::VertexIdx>   boundary;
        std::map<DCEL::VertexIdx, int> inside;
        uint32_t                       LEN = 0;
        for (size_t i = 0; i < output.Positions.size(); ++i) {
            DCEL::Vertex v = links.GetVertex(i);
            if (v.IsSide()) {
                boundary.push_back(i);
                inside[i]   = -1;
                glm::vec3 n = glm::normalize(output.Positions[i]);
                float     u = 0.5 * n.x + 0.5, v = 0.5 * n.y + 0.5;
                output.TexCoords.push_back(glm::vec2(u, v));
            } else {
                inside[i] = LEN++;
                output.TexCoords.push_back(glm::vec2(0.0, 0.0));
            }
        }
        // construct A,U2,V2
        float *A = new float[LEN * LEN], *U1 = new float[LEN], *U2 = new float[LEN], *V1 = new float[LEN], *V2 = new float[LEN];
        memset(A, 0, sizeof(float) * LEN * LEN);
        memset(U1, 0, sizeof(float) * LEN);
        memset(U2, 0, sizeof(float) * LEN);
        memset(V1, 0, sizeof(float) * LEN);
        memset(V2, 0, sizeof(float) * LEN);
        for (auto i : inside) {
            if (i.second >= 0) {
                DCEL::Vertex                 v         = links.GetVertex(i.first);
                std::vector<DCEL::VertexIdx> neighbors = v.GetNeighbors();
                size_t                       len       = neighbors.size();
                A[i.second * LEN + i.second]           = 1;
                for (size_t j = 0; j < len; ++j) {
                    if (inside[neighbors[j]] < 0) {
                        U2[i.second] += output.TexCoords[neighbors[j]].x * (float) 1 / len;
                        V2[i.second] += output.TexCoords[neighbors[j]].y * (float) 1 / len;
                    } else
                        A[i.second * LEN + inside[neighbors[j]]] = -(float) 1 / len;
                }
            }
        }
        // At=b
        for (uint32_t iter = 0; iter < numIterations; ++iter) {
            for (size_t y = 0; y < LEN; ++y) {
                U1[y] = U2[y];
                V1[y] = V2[y];
                for (size_t x = 0; x < LEN; ++x) {
                    if (A[y * LEN + x] == 1) continue;
                    if (A[y * LEN + x]) {
                        U1[y] -= A[y * LEN + x] * U1[x];
                        V1[y] -= A[y * LEN + x] * V1[x];
                    }
                }
            }
        }
        for (auto i : inside) {
            if (i.second >= 0)
                output.TexCoords[i.first] = glm::vec2(U1[i.second], V1[i.second]);
        }
        delete[] A;
        delete[] U1;
        delete[] U2;
        delete[] V1;
        delete[] V2;
        return;
    }

    /******************* 3. Mesh Simplification *****************/
    void SimplifyMesh(Engine::SurfaceMesh const & input, Engine::SurfaceMesh & output, float valid_pair_threshold, float simplification_ratio) {
        // your code here
        output = input;
        DCEL links;
        links.AddFaces(input.Indices);
        if (! links.IsValid()) {
            printf("3 : invalid");
            return;
        }
        std::vector<glm::mat4>                                q;
        std::set<std::pair<DCEL::VertexIdx, DCEL::VertexIdx>> pairs_set;
        std::vector<std::pair<DCEL::VertexIdx, DCEL::VertexIdx>> pairs;
        for (size_t i = 0; i < output.Positions.size(); ++i) {
            DCEL::Vertex vertex = links.GetVertex(i);
            // step 1
            std::vector<DCEL::Triangle const *> faces = vertex.GetFaces();
            glm::mat4                           Q(0);
            for (auto f : faces) {
                glm::vec3 p1 = output.Positions[f->Edges(0)->To()], p2 = output.Positions[f->Edges(1)->To()], p3 = output.Positions[f->Edges(2)->To()];
                float     a = (p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y);
                float     b = (p3.x - p1.x) * (p2.z - p1.z) - (p2.x - p1.x) * (p3.z - p1.z);
                float     c = (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y);
                float     l = sqrt(a * a + b * b + c * c);
                a /= l;
                b /= l;
                c /= l;
                float     d = -(a * p1.x + b * p1.y + c * p1.z);
                glm::vec4 p(a, b, c, d);
                glm::mat4 tmpQ(a * p, b * p, c * p, d * p);
                Q = Q + tmpQ;
            }
            q.push_back(Q);
            // step 2
            std::vector<DCEL::VertexIdx> neighbors = vertex.GetNeighbors();
            size_t                       len       = neighbors.size();
            for (size_t j = 0; j < len; ++j) {
                if (pairs_set.find(std::make_pair(neighbors[j], i)) != pairs_set.end()) continue; 
                pairs_set.insert(std::make_pair(i, neighbors[j]));
                pairs.push_back(std::make_pair(i, neighbors[j]));
            }
            for (size_t j = i + 1; j < output.Positions.size(); ++j) {
                if (i == j || pairs_set.find(std::make_pair(i, j)) != pairs_set.end()) continue;
                if (glm::distance(output.Positions[i], output.Positions[j]) < valid_pair_threshold) {
                    pairs_set.insert(std::make_pair(j, i));
                    pairs.push_back(std::make_pair(j, i));
                }
            }
        }
        // step 3&4
        std::vector<float> errors;
        for (auto p : pairs) {
            DCEL::VertexIdx v1 = p.first, v2 = p.second;
            glm::mat4       Q  = q[v1] + q[v2];
            glm::mat4       Q_ = Q;
            Q_[0][3]           = 0;
            Q_[1][3]           = 0;
            Q_[2][3]           = 0;
            Q_[3][3]           = 1;
            glm::vec4 v        = glm::inverse(Q_)[3];
            float     dis1 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[v1]), dis2 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[v2]);
            if (std::max(dis1, dis2) >= 2 * glm::distance(output.Positions[v1], output.Positions[v2]))
                v = glm::vec4(output.Positions[v1], 1);
            glm::vec4 tmp(glm::dot(v, Q[0]), glm::dot(v, Q[1]), glm::dot(v, Q[2]), glm::dot(v, Q[3]));
            float     error = glm::dot(v, tmp);
            errors.push_back(error);
        }
        // step 5
        size_t cnt = 0, loop_num = input.Positions.size() * (1 - simplification_ratio);
        while (cnt++ < loop_num) {
            float  MIN = 1000;
            size_t pos = 0;
            for (size_t i = 0; i < pairs.size(); ++i) {
                if (pairs[i].first == pairs[i].second) continue;
                if (errors[i] < MIN) {
                    MIN = errors[i];
                    pos = i;
                }
            }
            DCEL::VertexIdx v1 = pairs[pos].first, v2 = pairs[pos].second;
            pairs[pos].second = v1;
            errors[pos]       = 1000;
            glm::mat4 Q       = q[v1] + q[v2];
            q[v1]             = Q;
            glm::mat4 Q_      = Q;
            Q_[0][3]          = 0;
            Q_[1][3]          = 0;
            Q_[2][3]          = 0;
            Q_[3][3]          = 1;
            glm::vec4 v       = glm::inverse(Q_)[3];
            float     dis1 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[v1]), dis2 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[v2]);
            if (std::max(dis1, dis2) >= 2 * glm::distance(output.Positions[v1], output.Positions[v2]))
                v = glm::vec4(output.Positions[v1], 1);
            output.Positions[v1] = glm::vec3(v.x, v.y, v.z);
            for (size_t i = 0; i < pairs.size(); ++i) {
                if (pairs[i].first == pairs[i].second) continue;
                if (pairs[i].first == v2) {
                    pairs[i].first = v1;
                    glm::mat4 Q  = q[v1] + q[pairs[i].second];
                    glm::mat4 Q_ = Q;
                    Q_[0][3]     = 0;
                    Q_[1][3]     = 0;
                    Q_[2][3]     = 0;
                    Q_[3][3]     = 1;
                    glm::vec4 v  = glm::inverse(Q_)[3];
                    float     dis1 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[v1]), dis2 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[pairs[i].second]);
                    if (std::max(dis1, dis2) >= 2 * glm::distance(output.Positions[v1], output.Positions[pairs[i].second]))
                        v = glm::vec4(output.Positions[pairs[i].second], 1);
                    glm::vec4 tmp(glm::dot(v, Q[0]), glm::dot(v, Q[1]), glm::dot(v, Q[2]), glm::dot(v, Q[3]));
                    float     error = glm::dot(v, tmp);
                    errors[i]       = error;
                }
                if (pairs[i].second == v2) {
                    pairs[i].second = v1;
                    glm::mat4 Q  = q[v1] + q[pairs[i].first];
                    glm::mat4 Q_ = Q;
                    Q_[0][3]     = 0;
                    Q_[1][3]     = 0;
                    Q_[2][3]     = 0;
                    Q_[3][3]     = 1;
                    glm::vec4 v  = glm::inverse(Q_)[3];
                    float     dis1 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[v1]), dis2 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[pairs[i].first]);
                    if (std::max(dis1, dis2) >= 2 * glm::distance(output.Positions[v1], output.Positions[pairs[i].first]))
                        v = glm::vec4(output.Positions[pairs[i].first], 1);
                    glm::vec4 tmp(glm::dot(v, Q[0]), glm::dot(v, Q[1]), glm::dot(v, Q[2]), glm::dot(v, Q[3]));
                    float     error = glm::dot(v, tmp);
                    errors[i]       = error;
                }
                if (pairs[i].first == v1) {
                    glm::mat4 Q  = q[v1] + q[pairs[i].second];
                    glm::mat4 Q_ = Q;
                    Q_[0][3]     = 0;
                    Q_[1][3]     = 0;
                    Q_[2][3]     = 0;
                    Q_[3][3]     = 1;
                    glm::vec4 v  = glm::inverse(Q_)[3];
                    float     dis1 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[v1]), dis2 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[pairs[i].second]);
                    if (std::max(dis1, dis2) >= 2 * glm::distance(output.Positions[v1], output.Positions[pairs[i].second]))
                        v = glm::vec4(output.Positions[pairs[i].second], 1);
                    glm::vec4 tmp(glm::dot(v, Q[0]), glm::dot(v, Q[1]), glm::dot(v, Q[2]), glm::dot(v, Q[3]));
                    float     error = glm::dot(v, tmp);
                    errors[i]       = error;
                }
                if (pairs[i].second == v1) {
                    glm::mat4 Q  = q[v1] + q[pairs[i].first];
                    glm::mat4 Q_ = Q;
                    Q_[0][3]     = 0;
                    Q_[1][3]     = 0;
                    Q_[2][3]     = 0;
                    Q_[3][3]     = 1;
                    glm::vec4 v  = glm::inverse(Q_)[3];
                    float     dis1 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[v1]), dis2 = glm::distance(glm::vec3(v.x, v.y, v.z), output.Positions[pairs[i].first]);
                    if (std::max(dis1, dis2) >= 2 * glm::distance(output.Positions[v1], output.Positions[pairs[i].first]))
                        v = glm::vec4(output.Positions[pairs[i].first], 1);
                    glm::vec4 tmp(glm::dot(v, Q[0]), glm::dot(v, Q[1]), glm::dot(v, Q[2]), glm::dot(v, Q[3]));
                    float     error = glm::dot(v, tmp);
                    errors[i]       = error;
                }
            }
            for (size_t i = 0; i < output.Indices.size(); ++i) {
                if (output.Indices[i] == v2)
                    output.Indices[i] = v1;
            }
        }
    }

    /******************* 4. Mesh Smoothing *****************/
    void SmoothMesh(Engine::SurfaceMesh const & input, Engine::SurfaceMesh & output, std::uint32_t numIterations, float lambda, bool useUniformWeight) {
        // your code here
        output = input;
        DCEL links;
        links.AddFaces(output.Indices);
        if (! links.IsValid()) {
            printf("4 : invalid");
            return;
        }
        for (uint32_t n = 0; n < numIterations; ++n) {
            std::vector<glm::vec3> prepos = output.Positions;
            if (useUniformWeight) {
                for (std::size_t i = 0; i < output.Positions.size(); ++i) {
                    DCEL::Vertex          v         = links.GetVertex(i);
                    std::vector<uint32_t> neighbors = v.GetNeighbors();
                    int                   len       = neighbors.size();
                    glm::vec3             sum(0, 0, 0);
                    for (int j = 0; j < len; ++j)
                        sum = sum + prepos[neighbors[j]];
                    sum                 = sum / (float) len;
                    output.Positions[i] = (1 - lambda) * prepos[i] + lambda * sum;
                }
            } else {
                for (std::size_t i = 0; i < output.Positions.size(); ++i) {
                    DCEL::Vertex                        v         = links.GetVertex(i);
                    std::vector<uint32_t>               neighbors = v.GetNeighbors();
                    std::vector<const DCEL::Triangle *> faces     = v.GetFaces();
                    int                                 len       = neighbors.size();
                    glm::vec3                           sum(0, 0, 0);
                    float                               sum_w = 0;
                    for (int j = 0; j < len; ++j) {
                        float w = 0;
                        for (auto f : faces) {
                            int                    index = (f->IndiceOfVertex(i) + 2) % 3;
                            DCEL::HalfEdge const * e     = f->Edges(index);
                            if (e->To() != neighbors[j]) {
                                e = e->NextEdge();
                                if (e->To() != neighbors[j])
                                    continue;
                                e     = e->NextEdge();
                                index = (index + 2) % 3;
                            }
                            DCEL::VertexIdx v1   = e->OppositeVertex();
                            float           cot1 = glm::cot(glm::angle(prepos[v1] - prepos[i], prepos[v1] - prepos[neighbors[j]]));
                            w += cot1;
                            if (f->HasOppositeFace(index)) {
                                DCEL::VertexIdx v2   = e->PairOppositeVertex();
                                float           cot2 = glm::cot(glm::angle(prepos[v1] - prepos[i], prepos[v1] - prepos[neighbors[j]]));
                                w += cot2;
                            }
                            w = w > 0 ? w : -w;
                        }
                        sum_w += w;
                        sum = sum + w * prepos[neighbors[j]];
                    }
                    sum                 = sum / sum_w;
                    output.Positions[i] = (1 - lambda) * prepos[i] + lambda * sum;
                }
            }
        }
    }

    /******************* 5. Marching Cubes *****************/
    void MarchingCubes(Engine::SurfaceMesh & output, const std::function<float(const glm::vec3 &)> & sdf, const glm::vec3 & grid_min, const float dx, const int n) {
        // your code here
        size_t                 cnt = 0;
        std::vector<glm::vec3> unit;
        unit.push_back(glm::vec3(1, 0, 0));
        unit.push_back(glm::vec3(0, 1, 0));
        unit.push_back(glm::vec3(0, 0, 1));
        for (int x = 0; x < n; x++)
            for (int y = 0; y < n; y++)
                for (int z = 0; z < n; z++) {
                    glm::vec3 pos(x * dx, y * dx, z * dx);
                    pos       = pos + grid_min;
                    int index = 0;
                    for (int i = 0; i <= 7; ++i) {
                        glm::vec3 p(pos[0] + (i & 1) * dx, pos[1] + (i >> 1 & 1) * dx, pos[2] + (i >> 2 & 1) * dx);
                        if (sdf(p) < 0)
                            index |= (1 << i);
                    }
                    uint32_t e           = c_EdgeStateTable[index];
                    size_t   indices[12] = { 0 };
                    for (int j = 0; j < 12; ++j) {
                        if (e & (1 << j)) {
                            glm::vec3 p0 = pos + dx * (j & 1) * unit[((j >> 2) + 1) % 3] + dx * (j >> 1 & 1) * unit[((j >> 2) + 2) % 3];
                            glm::vec3 p1 = p0 + dx * unit[j >> 2];
                            glm::vec3 grad0(sdf(p0 + unit[0] * dx) - sdf(p0), sdf(p0 + unit[1] * dx) - sdf(p0), sdf(p0 + unit[2] * dx) - sdf(p0));
                            glm::vec3 grad1(sdf(p1 + unit[0] * dx) - sdf(p1), sdf(p1 + unit[1] * dx) - sdf(p1), sdf(p1 + unit[2] * dx) - sdf(p1));
                            grad0          = grad0 / dx;
                            grad1          = grad1 / dx;
                            float     rate = sdf(p0) / (sdf(p0) - sdf(p1));
                            glm::vec3 p    = p0 + dx * unit[j >> 2] * rate;
                            glm::vec3 grad = grad0 + (grad1 - grad0) * rate;
                            auto      tmp  = std::find(output.Positions.begin(), output.Positions.end(), p);
                            if (tmp != output.Positions.end())
                                indices[j] = tmp - output.Positions.begin();
                            else {
                                output.Positions.push_back(p);
                                output.Normals.push_back(grad);
                                indices[j] = cnt++;
                            }
                        }
                    }
                    std::array<int, 16> tri = c_EdgeOrdsTable[index];
                    int                 tmp = 0;
                    while (tri[tmp] != -1 && tmp < 16) {
                        output.Indices.push_back(indices[tri[tmp]]);
                        ++tmp;
                    }
                }
    }
} // namespace VCX::Labs::GeometryProcessing