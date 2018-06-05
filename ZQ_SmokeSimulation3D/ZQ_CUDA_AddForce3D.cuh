#ifndef _ZQ_CUDA_ADD_FORCE_3D_CUH_
#define _ZQ_CUDA_ADD_FORCE_3D_CUH_

namespace ZQ_CUDA_AddForce3D
{
	__global__
	void Compute_Vorticity_Vector_Scale_Kernel(float* vortVector, float* vortScale, const float* u, const float* v, const float* w, const float deltah, const int width, const int height, const int depth);

	__global__
	void Compute_Vorticity_Vector_Scale_Kernel(const int z, float* vortVector, float* vortScale, const float* u, const float* v, const float* w, const float deltah, const int width, const int height, const int depth);

	__global__
	void Compute_Gradient_of_VorticityScale_Kernel(float* gradVort, float* vortScale, const int width, const int height, const int depth);
	
	__global__
	void Compute_Gradient_of_VorticityScale_Kernel(const int z, float* gradVort, float* vortScale, const int width, const int height, const int depth);
	
	__global__
	void Compute_Force_Kernel(float* force, const float* gradVort, const float* vortVector, const float* temperature, const float confineCoeff, 
					const float buoyCoeff, const float deltah, const float Tamb, const int width, const int height, const int depth);

	__global__
	void Compute_Force_Kernel(const int z, float* force, const float* gradVort, const float* vortVector, const float* temperature, const float confineCoeff, 
					const float buoyCoeff, const float deltah, const float Tamb, const int width, const int height, const int depth);
	
	__global__
	void AddForce_u_Kernel(float* mac_u, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth);

	__global__
	void AddForce_u_Kernel(const int z, float* mac_u, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth);

	__global__
	void AddForce_v_Kernel(float* mac_v, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth);

	__global__
	void AddForce_v_Kernel(const int z, float* mac_v, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth);

	__global__
	void AddForce_w_Kernel(float* mac_w, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth);

	__global__
	void AddForce_w_Kernel(const int z, float* mac_w, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth);	
	
	/**********************************************************/

	void cu_Compute_Vorticity_Vector_Scale(float* vortVector, float* vortScale, float* u, float* v, float* w, const float deltah, const int width, const int height, const int depth);

	void cu_Compute_Vorticity_Vector_Scale2(float* vortVector, float* vortScale, float* u, float* v, float* w, const float deltah, const int width, const int height, const int depth);

	void cu_Compute_Gradient_of_VorticityScale(float* gradVort, float* vortScale, const int width, const int height, const int depth);

	void cu_Compute_Gradient_of_VorticityScale2(float* gradVort, float* vortScale, const int width, const int height, const int depth);

	void cu_Compute_Force(float* force, const float* gradVort, const float* vortVector, const float* temperature, const float confineCoeff, const float buoyCoeff, 
			const float deltah, const float Tamb, const int width, const int height, const int depth);

	void cu_Compute_Force2(float* force, const float* gradVort, const float* vortVector, const float* temperature, const float confineCoeff, const float buoyCoeff, 
			const float deltah, const float Tamb, const int width, const int height, const int depth);

	void cu_AddForce_u_v_w(float* mac_u, float* mac_v, float* mac_w, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth);

	void cu_AddForce_u_v_w2(float* mac_u, float* mac_v, float* mac_w, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth);

	void cu_AddForce3D(float* mac_u, float* mac_v, float* mac_w, const bool* occupy, const float* temperature, const float deltah, const float deltat, 
						const float buoyCoeff, const float confineCoeff, const float Tamb, const int width, const int height, const int depth);

	void cu_AddForce3D2(float* mac_u, float* mac_v, float* mac_w, const bool* occupy, const float* temperature, const float deltah, const float deltat, 
						const float buoyCoeff, const float confineCoeff, const float Tamb, const int width, const int height, const int depth);
		
}

#endif