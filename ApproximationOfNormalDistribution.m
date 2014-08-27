function N = ApproximationOfNormalDistribution(x)
	d1 = 0.0498673470;
	d2 = 0.0211410061;
	d3 = 0.0032776263;
	d4 = 0.0000380036;
	d5 = 0.0000488906;
	d6 = 0.0000053830;
    
    if x >= 0
        y = 1 + d1 * x + d2 * x^2 + d3 * x^3 + d4 * x^4 + d5 *x^5 + d6 * x^6;
        N = 1 - 0.5 * y^(-16);
    else
        x = -1 * x;
        n = ApproximationOfNormalDistribution(x);
        N = 1 - n;
    end
end