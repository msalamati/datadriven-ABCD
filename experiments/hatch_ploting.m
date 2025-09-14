clear all;close all
%parameters of the first patch object
obj(1).X = 1 + [0, 0 , 1 , 1 ];
obj(1).Y = 1 + [0, 1 , 1 , 0 ];
obj(1).FaceColor	= 'b';
obj(1).FaceAlpha	= 0.2;
HatchColor{1}		= [0.4,0.4,0.4];
HatchAngle{1}		= 20;
HatchType{1}		= 'single';
HatchDensity{1}		= 20;
%Parameters of the first patch object
obj(2).X = 0.5*[0, 0 , 1 , 1 ];
obj(2).Y = 0.5*[0, 1 , 1 , 0 ];
obj(2).FaceColor	= 'g';
obj(2).FaceAlpha	= 0.2;
HatchColor{2}		= [1,0,0];
HatchAngle{2}		= 20;
HatchDensity{2}		= 10;
HatchType{2}		= 'cross';
%make the patches
for ii=1:2
	h(ii) 		= patch(obj(ii));
end 
%hatch the patches
for ii=1:2
	hHatch(ii) 	= hatchfill2(h(ii),HatchType{ii},...
					'HatchAngle',HatchAngle{ii},...
					'HatchDensity',HatchDensity{ii},...
					'HatchColor',HatchColor{ii});
end
%Make the legend
Legend  			= {'blue square','green square'};
[legend_h,object_h,plot_h,text_str] = legendflex(h,Legend,...
 						'anchor',{'se','se'},...  
						'buffer',[-90, 60 ],... 
						'fontSize',20);
%object_h is the handle for the lines, patches and text objects
%hatch the legend patches to match the patches 
for ii=1:2
	hHatch(ii+2) 	= hatchfill2(object_h(ii+2),HatchType{ii}, ...
					'HatchAngle',HatchAngle{ii},...
					'HatchDensity',HatchDensity{ii},...
					'HatchColor',HatchColor{ii});
end
%I couldn't use hatchfill2 to set the FaceAlpha's so I did it manually
for ii= 1:2
	object_h(ii+2).FaceAlpha = obj(ii).FaceAlpha;  
end
