#include <math.h>
#include <stdlib.h>
#include "colormap.h"

void ColorMap::jetColorMap(unsigned char *rgb,float value,float min,float max)
{
	unsigned char c1=144;
	float max4=(max-min)/4;
	value-=min;
	if(value==HUGE_VAL)
	{rgb[0]=rgb[1]=rgb[2]=255;}
	else if(value<0)
	{rgb[0]=rgb[1]=rgb[2]=0;}
	else if(value<max4)
	{rgb[0]=0;rgb[1]=0;rgb[2]=c1+(unsigned char)((255-c1)*value/max4);}
	else if(value<2*max4)
	{rgb[0]=0;rgb[1]=(unsigned char)(255*(value-max4)/max4);rgb[2]=255;}
	else if(value<3*max4)
	{rgb[0]=(unsigned char)(255*(value-2*max4)/max4);rgb[1]=255;rgb[2]=255-rgb[0];}
	else if(value<max)
	{rgb[0]=255;rgb[1]=(unsigned char)(255-255*(value-3*max4)/max4);rgb[2]=0;}
	else {rgb[0]=255;rgb[1]=rgb[2]=0;}
}

void ColorMap::hotColorMap(unsigned char *rgb,float value,float min,float max)
{
	float max3=(max-min)/3;
	value-=min;
	if(value==HUGE_VAL)
	{rgb[0]=rgb[1]=rgb[2]=255;}
	else if(value<0)
	{rgb[0]=rgb[1]=rgb[2]=0;}
	else if(value<max3)
	{rgb[0]=(unsigned char)(255*value/max3);rgb[1]=0;rgb[2]=0;}
	else if(value<2*max3)
	{rgb[0]=255;rgb[1]=(unsigned char)(255*(value-max3)/max3);rgb[2]=0;}
	else if(value<max)
	{rgb[0]=255;rgb[1]=255;rgb[2]=(unsigned char)(255*(value-2*max3)/max3);}
	else {rgb[0]=rgb[1]=rgb[2]=255;}
}

void ColorMap::coldColorMap(unsigned char *rgb,float value,float min,float max)
{
	float max3=(max-min)/3;
	value-=min;
	if(value==HUGE_VAL)
	{rgb[0]=rgb[1]=rgb[2]=255;}
	else if(value<0)
	{rgb[0]=rgb[1]=rgb[2]=0;}
	else if(value<max3)
	{rgb[0]=0;rgb[1]=0;rgb[2]=(unsigned char)(255*value/max3);}
	else if(value<2*max3)
	{rgb[0]=0;rgb[1]=(unsigned char)(255*(value-max3)/max3);rgb[2]=255;}
	else if(value<max)
	{rgb[0]=(unsigned char)(255*(value-2*max3)/max3);rgb[1]=255;rgb[2]=255;}
	else {rgb[0]=rgb[1]=rgb[2]=255;}
}

void ColorMap::blueColorMap(unsigned char *rgb,float value,float min,float max)
{
	value-=min;
	if(value==HUGE_VAL)
	{rgb[0]=rgb[1]=rgb[2]=255;}
	else if(value<0)
	{rgb[0]=rgb[1]=rgb[2]=0;}
	else if(value<max)
	{rgb[0]=0;rgb[1]=0;rgb[2]=(unsigned char)(255*value/max);}
	else {rgb[0]=rgb[1]=0;rgb[2]=255;}
}

void positiveColorMap(unsigned char *rgb,float value,float min,float max)

{
	value-=min;
	max-=min;
	value/=max;

	if(value<0){
		rgb[0]=rgb[1]=rgb[2]=0;
		return;
	}
	if(value>1){
		rgb[0]=rgb[1]=rgb[2]=255;
		return;
	}

	rgb[0]=192;rgb[1]=0;rgb[2]=0;
	rgb[0]+=(unsigned char)(63*value);
	rgb[1]+=(unsigned char)(255*value);
	if(value>0.5)
		rgb[2]+=(unsigned char)(255*2*(value-0.5));
}

void negativeColorMap(unsigned char *rgb,float value,float min,float max)
{
	value-=min;

	max-=min;

	value/=max;

	rgb[0]=0;rgb[1]=0;rgb[2]=0;

	if(value<0) return;

	if(value>1){
		rgb[1]=rgb[2]=255;
		return;
	}

	rgb[1]+=(unsigned char)(255*value);
	
	if(value > 0.5)
		rgb[2]+=(unsigned char)(255*2*(value-0.5));
}

void ColorMap::colorMap(unsigned char *rgb,float value,float min,float max)
{
	if(value>0) 
		positiveColorMap(rgb,value,0,max);
	else 
		negativeColorMap(rgb,value,min,0);
	/*
	if(value>0) 
		hotColorMap(rgb,value,min,max);
	else 
		coldColorMap(rgb,value,min,max);
	*/
}

void ColorMap::cyclicColorMap(unsigned char *rgb,float value,float min,float max)
{
	float max3=(max-min)/3;
	value-=(max-min)*(float)floor((value-min)/(max-min));
	if(value<max3)
	{rgb[0]=(unsigned char)(255-255*value/max3);rgb[1]=0;rgb[2]=255-rgb[0];}
	else if(value<2*max3)
	{rgb[0]=0;rgb[1]=(unsigned char)(255*(value-max3)/max3);rgb[2]=255-rgb[1];}
	else if(value<max)
	{rgb[0]=(unsigned char)(255*(value-2*max3)/max3);rgb[1]=255-rgb[0];rgb[2]=0;}

}
void ColorMap::randColorMap(unsigned char *rgb,float value,float min,float max)
{
	srand((int)(65000*(value-min)/(max-min)));
	rgb[0]=(unsigned char)(255*rand());
	rgb[1]=(unsigned char)(255*rand());
	rgb[2]=(unsigned char)(255*rand());
}

void ColorMap::grayColorMap(unsigned char *rgb,float value,float min,float max)
{
	max-=min;
	value-=min;
	rgb[0]=rgb[1]=rgb[2]=(unsigned char)(255*value/max);
}

colorMapFunc ColorMap::selectColorMap(int cmp)
{
	int max=7;
	cmp=abs(cmp)%max;
	switch(cmp){
	case 0:
		return colorMap;
		break;
	case 1:
		return hotColorMap;
		break;
	case 2:
		return coldColorMap;
		break;
	case 3:
		return jetColorMap;
		break;
	case 4:
		return cyclicColorMap;
		break;
	case 5:
		return grayColorMap;
		break;
	case 6:
		return blueColorMap;
		break;
	default:
		return colorMap;
	};
}
