package plugins.pikeja.colocalizationstudiomod;

import org.apache.poi.ss.formula.OperationEvaluationContext;
import org.apache.poi.ss.formula.eval.EvaluationException;
import org.apache.poi.ss.formula.eval.NumberEval;
import org.apache.poi.ss.formula.eval.OperandResolver;
import org.apache.poi.ss.formula.eval.ValueEval;
import org.apache.poi.ss.formula.functions.FreeRefFunction;

import flanagan.analysis.Stat;

public class Quantile implements FreeRefFunction {
	@Override
	public ValueEval evaluate(ValueEval[] args, OperationEvaluationContext ec) {		
		try {
			final ValueEval alpha = OperandResolver.getSingleValue(args[0], ec.getRowIndex(), ec.getColumnIndex()) ;
			final ValueEval k = OperandResolver.getSingleValue(args[1], ec.getRowIndex(), ec.getColumnIndex()) ;
			final ValueEval n = OperandResolver.getSingleValue(args[2], ec.getRowIndex(), ec.getColumnIndex()) ;			
			return new NumberEval(Quantile.compute(OperandResolver.coerceValueToDouble(alpha), OperandResolver.coerceValueToInt(k), OperandResolver.coerceValueToInt(n)));
		} catch (EvaluationException e) {
			return e.getErrorEval();
		}			 		
	}
	
	/**
	 * Original C++ implementation found at http://www.wilmott.com/messageview.cfm?catid=10&threadid=38771
	 * C# implementation found at http://weblogs.asp.net/esanchez/archive/2010/07/29/a-quick-and-dirty-implementation-of-excel-norminv-function-in-c.aspx
	 *    Compute the quantile function for the normal distribution.
	 *
	 *    For small to moderate probabilities, algorithm referenced
	 *    below is used to obtain an initial approximation which is
	 *    polished with a final Newton step.
	 *    For very large arguments, an algorithm of Wichura is used.
	 *
	 *  REFERENCE
	 *
	 *    Beasley, J. D. and S. G. Springer (1977).
	 *    Algorithm AS 111: The percentage points of the normal distribution,
	 *    Applied Statistics, 26, 118-121.
	 *
	 *     Wichura, M.J. (1988).
	 *     Algorithm AS 241: The Percentage Points of the Normal Distribution.
	 *     Applied Statistics, 37, 477-484.
	 * @param p
	 * @param mu
	 * @param sigma
	 * @return
	 */
	public static double compute(double alpha, int k, int n) {
		double n_1,quantile,pas,p_q,proba;		
		if(alpha < 0 || alpha > 1) 
			throw new RuntimeException("alpha must be bigger than 0 and smaller than 1");		
		if(k < 1 || k > n)
			throw new RuntimeException("k must be bigger than 1 and smaller than n");
		if(n < 1)
			throw new RuntimeException("n must be bigger than 1");
		if(alpha == 0)
			return Double.NEGATIVE_INFINITY;		
		if(alpha == 1)
			return Double.POSITIVE_INFINITY;		
		if(k == 1)//retour quantile sup
			{n_1=(double) 1/(n);
			quantile=NormInv.compute(Math.pow(alpha,n_1),0,1);
			return quantile;}
		if(k == n)//retour quantile inf
			{n_1=(double) 1/(n);
			quantile=NormInv.compute(1-Math.pow(1-alpha,n_1),0,1);
			return quantile;
			}
		//situation intermediaire on fait une boucle 
		pas=0.0001;quantile=-5;
		proba=1;double bin=0;
		while (proba>(1-alpha))
		{proba=0;p_q=Stat.normalCDF(0, 1, quantile);
			for (int j=0;j<n-k+1;j++)
			{
				bin=Stat.binomialCoeff(n, j);
				//bin= Binomial.compute (j,n);
				proba+=bin*Math.pow(1-Stat.normalCDF(0, 1, quantile), n-j)*Math.pow(Stat.normalCDF(0, 1, quantile), j);
			}					
			quantile+=pas;
		}		
		return quantile;		
	}
	 
}