function hp = TrackStatsPlotCI( x ,y, color, SmoothingWidth , Statisticly_Dependent_Coeff)
% used to plot moving-window-binned over x
% + coturn handle to plot (for legend etc) (without STD bars)

	if ~exist( 'SmoothingWidth', 'var' )
		SmoothingWidth = 0.5;
    end
    
    if ~exist( 'Statisticly_Dependent_Coeff', 'var' )
		Statisticly_Dependent_Coeff = 1;
    end
    
	if ~exist( 'color', 'var' )
		color = 'g';
    end

%     if false
% 	[xDataOut, yDataOut, yCI] = SteffenSmoothCI( x, y, SmoothingWidth );
% 	hp = plot( xDataOut(yCI<1 & yCI>0), yDataOut(yCI<1 & yCI>0), 'linewidth', 2, 'color', color, 'LineWidth', 2 );
% 	hold on
% 	ciplot_steffen( xDataOut(yCI<1 & yCI>0), yDataOut(yCI<1 & yCI>0)-yCI(yCI<1 & yCI>0)*Statisticly_Dependent_Coeff^.5, yDataOut(yCI<1 & yCI>0)+yCI(yCI<1 & yCI>0)*Statisticly_Dependent_Coeff^.5, color );
%     end
    
%% run this for plotting 
    if false %define upper and lower values
     SmoothingWidth = 0.1;
        high_bound=5.7;
        low_bound=4.8;
        x=Track_mean_nuc_shape; y=motility_measure2';
    [xDataOut, yDataOut, yCI] = SteffenSmoothCI( x, y, SmoothingWidth );
    hp = plot( xDataOut(xDataOut(:)<high_bound & xDataOut(:)>low_bound), yDataOut(xDataOut(:)<high_bound & xDataOut(:)>low_bound), 'linewidth', 2, 'color', color, 'LineWidth', 2 );
	hold on
	ciplot_steffen( xDataOut(xDataOut(:)<high_bound & xDataOut(:)>low_bound), yDataOut(xDataOut(:)<high_bound & xDataOut(:)>low_bound)-yCI(xDataOut(:)<high_bound & xDataOut(:)>low_bound)*Statisticly_Dependent_Coeff^.5, yDataOut(xDataOut(:)<high_bound & xDataOut(:)>low_bound)+yCI(xDataOut(:)<high_bound & xDataOut(:)>low_bound)*Statisticly_Dependent_Coeff^.5, color ); 

    end 
%%
%     [xDataOut, yDataOut, yCI] = SteffenSmoothCI( x, y, SmoothingWidth );
%     hp = plot( xDataOut(:), yDataOut(:), 'linewidth', 2, 'color', color, 'LineWidth', 2 );
% 	hold on
% 	ciplot_steffen( xDataOut(:), yDataOut(:)-yCI(:)*Statisticly_Dependent_Coeff^.5, yDataOut(:)+yCI(:)*Statisticly_Dependent_Coeff^.5, color );  
    
end

