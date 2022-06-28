function [ar,eG,el_shape,dG,axis,eV]=get_axis_of_shape_in_box(shapeinbox,ratio)

    [X,Y,Z]=meshgrid(1:size(shapeinbox,2),1:size(shapeinbox,1),1:size(shapeinbox,3));
    
    % find center of mass
	mass = sum(sum(sum( shapeinbox )));
    cx = sum(sum(sum( X.*shapeinbox ))) / mass;
    cy = sum(sum(sum( Y.*shapeinbox ))) / mass;
    cz = sum(sum(sum( Z.*shapeinbox ))) / mass;
        
    % correct coordinate system
	Xc		= (X - cx)*ratio(1);
	Yc		= (Y - cy)*ratio(2);
    Zc		= (Z - cz)*ratio(3);
	mass	= mass * prod(ratio); % auch korrigieren by volume element!
    
    % second area moments
	Ixx = sum(sum(sum( shapeinbox.*Xc.*Xc )));
	Ixy = sum(sum(sum( shapeinbox.*Xc.*Yc )));
	Iyy = sum(sum(sum( shapeinbox.*Yc.*Yc )));
    Ixz = sum(sum(sum( shapeinbox.*Xc.*Zc )));
    Iyz = sum(sum(sum( shapeinbox.*Yc.*Zc )));
    Izz = sum(sum(sum( shapeinbox.*Zc.*Zc )));

	% moment matrix
	dG = [Ixx, Ixy, Ixz; Ixy, Iyy, Iyz; Ixz, Iyz, Izz] / mass;%* prod(ratio); % *prod ratio weils ansonsten zweimal rausnormiert ist
    eG = eig(dG);
    [eV,~] = eig(dG);
    
    ar=(max(eG)/min(eG))^.5;
    a = (5*eG(1))^0.5;
    b = (5*eG(2))^0.5;
    c = (5*eG(3))^0.5;
	
    axis=[a,b,c];
    el_shape = 4*pi* ( ((a*b)^(8/5)+(a*c)^(8/5)+(c*b)^(8/5))/3 )^(5/8)  / (4/3*pi*a*b*c)^(2/3);
	
	
end