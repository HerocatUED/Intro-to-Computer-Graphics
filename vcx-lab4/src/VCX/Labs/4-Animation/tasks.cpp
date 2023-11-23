#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <spdlog/spdlog.h>
#include <iostream>
#include <cmath>
#include "Labs/4-Animation/tasks.h"
#include "IKSystem.h"
#include "CustomFunc.inl"


namespace VCX::Labs::Animation {
    void ForwardKinematics(IKSystem & ik, int StartIndex) {
        if (StartIndex == 0) {
            ik.JointGlobalRotation[0] = ik.JointLocalRotation[0];
            ik.JointGlobalPosition[0] = ik.JointLocalOffset[0];
            StartIndex                = 1;
        }

        for (int i = StartIndex; i < ik.JointLocalOffset.size(); i++) {
            // your code here: forward kinematics
            ik.JointGlobalRotation[i] = glm::normalize(ik.JointGlobalRotation[i - 1] * ik.JointLocalRotation[i]);
            ik.JointGlobalPosition[i] = ik.JointGlobalPosition[i - 1] + glm::rotate(ik.JointGlobalRotation[i - 1], ik.JointLocalOffset[i]);
        }
    }
         
    void InverseKinematicsCCD(IKSystem & ik, const glm::vec3 & EndPosition, int maxCCDIKIteration, float eps) {
        ForwardKinematics(ik, 0);
        // These functions will be useful: glm::normalize, glm::rotation, glm::quat * glm::quat
        for (int CCDIKIteration = 0; CCDIKIteration < maxCCDIKIteration && glm::l2Norm(ik.EndEffectorPosition() - EndPosition) > eps; CCDIKIteration++) {
            // your code here: ccd ik
            int end = ik.NumJoints() - 1;
            for (int i = end - 1; i >= 0; --i) {
                glm::vec3 toEnd          = ik.JointGlobalPosition[end] - ik.JointGlobalPosition[i];
                glm::vec3 toTarget       = EndPosition - ik.JointGlobalPosition[i];
                glm::quat Rotation       = glm::rotation(glm::normalize(toEnd), glm::normalize(toTarget));
                ik.JointLocalRotation[i] = glm::normalize(Rotation * ik.JointLocalRotation[i]);
                ForwardKinematics(ik, i);
            }
        }
    }

    void InverseKinematicsFABR(IKSystem & ik, const glm::vec3 & EndPosition, int maxFABRIKIteration, float eps) {
        ForwardKinematics(ik, 0);
        int                    nJoints = ik.NumJoints();
        std::vector<glm::vec3> backward_positions(nJoints, glm::vec3(0, 0, 0)), forward_positions(nJoints, glm::vec3(0, 0, 0));
        for (int IKIteration = 0; IKIteration < maxFABRIKIteration && glm::l2Norm(ik.EndEffectorPosition() - EndPosition) > eps; IKIteration++) {
            // task: fabr ik
            // backward update
            glm::vec3 next_position         = EndPosition;
            backward_positions[nJoints - 1] = EndPosition;

            for (int i = nJoints - 2; i >= 0; i--) {
                // your code here
                glm::vec3 direction = glm::normalize(ik.JointGlobalPosition[i] - next_position);
                next_position += direction * ik.JointOffsetLength[i + 1];
                backward_positions[i] = next_position;
            }

            // forward update
            glm::vec3 now_position = ik.JointGlobalPosition[0];
            forward_positions[0]   = ik.JointGlobalPosition[0];
            for (int i = 0; i < nJoints - 1; i++) {
                // your code here
                glm::vec3 direction = glm::normalize(backward_positions[i + 1] - now_position);
                now_position += direction * ik.JointOffsetLength[i + 1];
                forward_positions[i + 1] = now_position;
            }
            ik.JointGlobalPosition = forward_positions; // copy forward positions to joint_positions
        }

        // Compute joint rotation by position here.
        for (int i = 0; i < nJoints - 1; i++) {
            ik.JointGlobalRotation[i] = glm::rotation(glm::normalize(ik.JointLocalOffset[i + 1]), glm::normalize(ik.JointGlobalPosition[i + 1] - ik.JointGlobalPosition[i]));
        }
        ik.JointLocalRotation[0] = ik.JointGlobalRotation[0];
        for (int i = 1; i < nJoints - 1; i++) {
            ik.JointLocalRotation[i] = glm::inverse(ik.JointGlobalRotation[i - 1]) * ik.JointGlobalRotation[i];
        }
        ForwardKinematics(ik, 0);
    }

    IKSystem::Vec3ArrPtr IKSystem::BuildCustomTargetPosition() {
        // get function from https://www.wolframalpha.com/input/?i=Albert+Einstein+curve
        int nums      = 500;
        //int nums      = 5000;
        using Vec3Arr = std::vector<glm::vec3>;
        std::shared_ptr<Vec3Arr> custom(new Vec3Arr(nums));
        int                      index = 0;
        float                    a = 0.1, pi = acos(-1);
        for (int i = 0; i < nums; i++) {
            // float x_val = 1.5e-3f * custom_x(92 * glm::pi<float>() * i / nums);
            // float y_val = 1.5e-3f * custom_y(92 * glm::pi<float>() * i / nums);
            // if (std::abs(x_val) < 1e-3 || std::abs(y_val) < 1e-3) continue;
            //(*custom)[index++] = glm::vec3(1.6f - x_val, 0.0f, y_val - 0.2f);

            float x_val        = a * (6 * cos(2 * pi * i / nums) - cos(6 * 2 * pi * i / nums));
            float y_val        = a * (6 * sin(2 * pi * i / nums) - sin(6 * 2 * pi * i / nums));
            (*custom)[index++] = glm::vec3(x_val, 0.0f, y_val);
        }
        custom->resize(index);
        return custom;
    }

    void AdvanceMassSpringSystem(MassSpringSystem & system, float const dt) {
        // your code here: rewrite following code
        int const   steps = 50;
        float const ddt   = dt / steps;
        int         cnt   = system.Positions.size();
        for (std::size_t s = 0; s < steps; s++) {
            //  float                  E = 0;
            //  float                  inertia = 0;
            std::vector<glm::vec3> forces(system.Positions.size(), glm::vec3(0));
            std::vector<glm::vec3> XY(cnt, glm::vec3(0)); // X - Y
            for (auto const spring : system.Springs) {
                auto const      p0  = spring.AdjIdx.first;
                auto const      p1  = spring.AdjIdx.second;
                glm::vec3 const x01 = system.Positions[p1] - system.Positions[p0];
                glm::vec3 const v01 = system.Velocities[p1] - system.Velocities[p0];
                glm::vec3 const e01 = glm::normalize(x01);
                glm::vec3       f   = (system.Stiffness * (glm::length(x01) - spring.RestLength) + system.Damping * glm::dot(v01, e01)) * e01;
                forces[p0] += f;
                forces[p1] -= f;
                // float delta_x = glm::length(x01) - spring.RestLength;
                // float energe  = system.Stiffness * delta_x * delta_x / 2;
                //  E += energe;
            }
            for (int i = 0; i < cnt; i++) {
                XY[i] = -ddt * (system.Velocities[i] + ddt * glm::vec3(0, -system.Gravity, 0)) / system.Mass;
                // inertia += glm::length(X[i] - Y[i]) / (2 * ddt * ddt);
            }
            // float gx = inertia + E;

            Eigen::MatrixXf grad_g(3 * cnt, 1);
            float           hm = system.Mass / (ddt * ddt);
            for (int i = 0; i < cnt; i++) {
                grad_g(3 * i, 0)     = hm * XY[i].x - forces[i].x;
                grad_g(3 * i + 1, 0) = hm * XY[i].y - forces[i].y;
                grad_g(3 * i + 2, 0) = hm * XY[i].z - forces[i].z;
            }

            Eigen::SparseMatrix<float> Hessian(3 * cnt, 3 * cnt);
            for (int i = 0; i < 3 * cnt; ++i)
                Hessian.coeffRef(i, i) = hm;
            for (auto const spring : system.Springs) {
                auto const      p0      = spring.AdjIdx.first;
                auto const      p1      = spring.AdjIdx.second;
                glm::vec3 const x01     = system.Positions[p1] - system.Positions[p0];
                float           len     = glm::length(system.Positions[p1] - system.Positions[p0]);
                float           inv_len = 1.0f / len;
                float           xy      = system.Stiffness * x01.x * x01.y * inv_len * inv_len;
                float           xz      = system.Stiffness * x01.x * x01.z * inv_len * inv_len;
                float           yz      = system.Stiffness * x01.y * x01.z * inv_len * inv_len;
                float           xx      = system.Stiffness * inv_len * (len - spring.RestLength + x01.x * x01.x * inv_len);
                float           yy      = system.Stiffness * inv_len * (len - spring.RestLength + x01.y * x01.y * inv_len);
                float           zz      = system.Stiffness * inv_len * (len - spring.RestLength + x01.z * x01.z * inv_len);

                Hessian.coeffRef(3 * p0, 3 * p0) += xx;
                Hessian.coeffRef(3 * p0, 3 * p0 + 1) += xy;
                Hessian.coeffRef(3 * p0, 3 * p0 + 2) += xz;
                Hessian.coeffRef(3 * p0 + 1, 3 * p0 + 1) += yy;
                Hessian.coeffRef(3 * p0 + 1, 3 * p0) += xy;
                Hessian.coeffRef(3 * p0 + 1, 3 * p0 + 2) += yz;
                Hessian.coeffRef(3 * p0 + 2, 3 * p0 + 2) += zz;
                Hessian.coeffRef(3 * p0 + 2, 3 * p0) += xz;
                Hessian.coeffRef(3 * p0 + 2, 3 * p0 + 1) += yz;

                Hessian.coeffRef(3 * p1, 3 * p1) += xx;
                Hessian.coeffRef(3 * p1, 3 * p1 + 1) += xy;
                Hessian.coeffRef(3 * p1, 3 * p1 + 2) += xz;
                Hessian.coeffRef(3 * p1 + 1, 3 * p1 + 1) += yy;
                Hessian.coeffRef(3 * p1 + 1, 3 * p1) += xy;
                Hessian.coeffRef(3 * p1 + 1, 3 * p1 + 2) += yz;
                Hessian.coeffRef(3 * p1 + 2, 3 * p1 + 2) += zz;
                Hessian.coeffRef(3 * p1 + 2, 3 * p1) += xz;
                Hessian.coeffRef(3 * p1 + 2, 3 * p1 + 1) += yz;
            }
            // x
            auto            solver = Eigen::SimplicialLLT<Eigen::SparseMatrix<float>>(Hessian);
            Eigen::MatrixXf delta_x(3 * cnt, 1);
            delta_x = solver.solve(grad_g);
            for (int i = 0; i < cnt; ++i) {
                if (system.Fixed[i]) continue;
                system.Positions[i] -= glm::vec3(delta_x(3 * i, 0), delta_x(3 * i + 1, 0), delta_x(3 * i + 2, 0));
            }
            // v
            for (auto const spring : system.Springs) {
                auto const      p0  = spring.AdjIdx.first;
                auto const      p1  = spring.AdjIdx.second;
                glm::vec3 const x01 = system.Positions[p1] - system.Positions[p0];
                glm::vec3 const v01 = system.Velocities[p1] - system.Velocities[p0];
                glm::vec3 const e01 = glm::normalize(x01);
                glm::vec3       f   = (system.Stiffness * (glm::length(x01) - spring.RestLength) + system.Damping * glm::dot(v01, e01)) * e01;
                forces[p0] += f;
                forces[p1] -= f;
            }
            for (int i = 0; i < cnt; ++i) {
                if (system.Fixed[i]) continue;
                system.Velocities[i] += (glm::vec3(0, -system.Gravity, 0) + forces[i] / system.Mass) * ddt;
            }
        }
    }
}
