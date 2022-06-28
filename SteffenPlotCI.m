function [hp, xDataOut, yDataOut, yCI] = SteffenPlotCI( x ,y, color, SmoothingWidth )
% used to plot moving-window-binned over x
% + confidence interval
% return handle to plot (for legend etc) (without STD bars)

	if ~exist( 'SmoothingWidth', 'var' )
		SmoothingWidth = 0.2;
	end

	[xDataOut, yDataOut, yCI] = SteffenSmoothCI( x, y, SmoothingWidth );
	hp = plot( xDataOut, yDataOut, 'linewidth', 2, 'color', color, 'LineWidth', 2 )
	hold on
	ciplot_steffen( xDataOut, yDataOut-yCI, yDataOut+yCI, color );

end

