function [xDataOut, yDataOut, yCI] = SteffenSmoothCI( xData, yData, WindowWidth, SampleNumber )

	% WindowWidth:
	% Objekte in +WindowWidth/2 bzw -WindowWidth/2 Entfernung 
	% werden mit 0.5 gewichtet


	% SampleNumber : how many output sample points should there be (equally
	% distributed bet xMin and xMax)
	% -> defaults to 100


	% check input
	if length(xData) ~= length(yData)
		error('SteffenSmooth: x and y input vectors do not have the same length!');
	end
	
	if ~exist('SampleNumber', 'var')
		SampleNumber = 100;
	end
	
	SampleNumber = round(SampleNumber);
	if SampleNumber < 1
		error('SteffenSmooth: Samplenumber must > 1...');
	end
	
	% kick out all the data points where (x or y) is NaN!
	keep = find( ~ (isnan(xData) | isnan(yData)) );
	xData = xData(keep);
	yData = yData(keep);
	

	% create the x output vector equally spaced from MinX to MaxX
	xDataOut	= min(xData) + (0:SampleNumber-1)/(SampleNumber-1)*(max(xData)-min(xData));

	% create the y output vectors
	yDataOut	= zeros( size(xDataOut) );
	yStd		= zeros( size(xDataOut) );
	yCI			= zeros( size(xDataOut) );
	
	% for confidence interval
	alpha		= 0.05;
	
	for j=1:SampleNumber
		% each data point weighed according to its distance to the current
		% output sample point
		InputWeights	= WeightWindow( xData-xDataOut(j) );
		TotalWeight		= sum( InputWeights );
		
		% weigh the y data
		yDataOut(j) = sum(yData.*InputWeights) / TotalWeight;
		
		% construct sth like the STD deviation of data around this sample
		% point - of course also weighted ...
		yStd(j)		= sqrt( sum( (yData-yDataOut(j)) .^ 2 .* InputWeights) / TotalWeight );
		% This gives narrower "std" for more values
		SEM			= yStd(j) / sqrt( TotalWeight );
		
		% quantil for the given alpha value (e.g. 0.05)
		% for large enough sample numbers, this approaches a value of ~2
		% (strictly spealking this only works for gaussian disrtb°s!)
		quantil		= tinv( 1-alpha/2, TotalWeight-1 );

		yCI(j)		= SEM * quantil;
	end
	
	% die weight funktion
	function w = WeightWindow( dx )
		% Objekte in +WindowWidth/2 bzw -WindowWidth/2 werden mit 0.5 gewichtet:		
		Sigma = WindowWidth/2.4;
		w = exp( - dx.^2 / (2*Sigma^2) );
	end

end