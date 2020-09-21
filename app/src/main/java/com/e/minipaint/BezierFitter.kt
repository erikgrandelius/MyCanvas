package com.e.minipaint

val Boolean.double
    get() = if (this) 1.0 else 0.0

class BezierFitter constructor(sampledX: Array<Number>, sampledY: Array<Number>, tol: Number = 1e-4){
    fun solveTridiagonal(A : List<DoubleArray>,  b : DoubleArray ): DoubleArray {
        /*
        * Implements the Thomas algorithm for solving the tri-diagonal linear system.
        * https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm.
        * Only stable in the diagonally dominant or positive definite case.
        *Ax = b, where A = diag(d_lower , d_diag, d_upper).
        *
        * Returns: x: Array<Double>
                   solution
        */
        println(A)

        val d_lower = A[0]
        val d_diag = A[1]
        val d_upper = A[2]
        val n = d_diag.size
        val d_upper_temp = DoubleArray(n-1)
        val b_temp = DoubleArray(n)

        d_upper_temp[0] = d_upper[0]/d_diag[0]
        b_temp[0] = b[0]/d_diag[0]
        for (i in 1 .. n - 2) {
            d_upper_temp[i] = d_upper[i]/(d_diag[i]-d_lower[i]*d_upper_temp[i-1])
            b_temp[i] = (b[i] - b_temp[i-1] * d_lower[i] )/(d_diag[i] - d_lower[i] * d_upper_temp[i-1])
            println(i.toString(i))}
        b[n-1] = (b[n-1] - b_temp[n-2] * d_lower[n-2])/(d_diag[n-1] - d_lower[n-1] * d_upper_temp[n-2])
        val x = DoubleArray(n)
        x[n-1] = b_temp[n-1]
        for (i in n-2 downTo 0) {
            x[i] = (b_temp[i] - d_upper[i] * x[i+1])/d_diag[i]
        }
        return x
    }

    fun getBezierCoef(x: DoubleArray, y: DoubleArray): List<DoubleArray>{
        val n = x.size - 1
        val C = listOf<DoubleArray>(DoubleArray(n-1, { i -> 1 + (i == n-1).double}),
            DoubleArray(n, {i -> 4 - 2 * (i == 0).double + 3 * (i == n).double}),
            DoubleArray(n- 1, {i -> 1.0})
        )
        val P_x = DoubleArray(n, { i -> 2 * (2 * x[i] + x[i + 1])})
        P_x[0] = x[0] + 2 * x[1]
        P_x[n-1] = 8 * x[n - 1] + x[n]
        val P_y = DoubleArray(n, {i -> 2 * (2 * y[i] + y[i + 1])})
        P_y[0] = y[0] + 2 * y[1]
        P_y[n-1] = 8 * y[n - 1] + y[n]
        val A_x = this.solveTridiagonal(C, P_x)
        val A_y = this.solveTridiagonal(C, P_y)
        val B_x = DoubleArray(n, { i -> if (i < n-1) 2 * x[i + 1] - A_x[i + 1] else (A_x[n - 1] + x[n]) / 2})
        val B_y = DoubleArray(n, { i -> if (i < n-1) 2 * y[i + 1] - A_y[i + 1] else (A_y[n - 1] + y[n]) / 2})

        return listOf(A_x, A_y, B_x, B_y)
    }
}