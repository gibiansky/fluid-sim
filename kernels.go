package main

import "math"

/** Smoothing Kernels:
  - General Kernel: used for most SPH computations
  - Pressure Kernel: used for computing pressure
  - Viscosity Kernel: used for computing acceleration due to viscosity
  **/

// General kernel
func smoothingKernel(radius float64) float64 {
	return (315 / (64 * math.Pi * math.Pow(kernelRadius, 9))) * math.Pow(math.Pow(kernelRadius, 2)-math.Pow(radius, 2), 3)
}

// Derivative of general kernel
func derivSmoothingKernel(radius float64) float64 {
	return -(315 * 3 / (64 * math.Pi * math.Pow(kernelRadius, 9))) * math.Pow(math.Pow(kernelRadius, 2)-math.Pow(radius, 2), 2) * 2 * radius
}

// Laplacian of general smoothing kernel
func laplaceSmoothingKernel(radius float64) float64 {
	return -945 / 32 / math.Pow(kernelRadius, 9) / math.Pi * (3*math.Pow(kernelRadius, 2) - 7*math.Pow(radius, 2)) * (math.Pow(kernelRadius, 2) - math.Pow(radius, 2))
}

// Pressure kernel
func pressureKernel(radius float64) float64 {
	return 15 / (math.Pi * math.Pow(kernelRadius, 6)) * math.Pow(kernelRadius-radius, 3)
}

// Gradient of pressure kernel
func derivPressureKernel(radius float64) float64 {
	return -45 / (math.Pi * math.Pow(kernelRadius, 6)) * math.Pow(kernelRadius-radius, 2)
}

// Viscosity kernel
func viscosityKernel(radius float64) float64 {
	ker := -(math.Pow(radius, 3) / (2 * math.Pow(kernelRadius, 3))) + (math.Pow(radius, 2) / (math.Pow(kernelRadius, 2))) + (kernelRadius / (2 * radius)) - 1
	return ker * 15 / 2 / math.Pow(kernelRadius, 3) / math.Pi
}

// Laplacian of viscosity kernel
func laplaceViscosityKernel(radius float64) float64 {
	ker := 6 * (kernelRadius - radius) / math.Pow(kernelRadius, 3)
	return ker * 15 / 2 / math.Pow(kernelRadius, 3) / math.Pi
}
