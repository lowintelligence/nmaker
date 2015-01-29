/*
 * =====================================================================================
 *
 *       Filename:  ppverify.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/25/2015 07:48:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Cao Zongyan (), zycao@sccas.cn
 *        Company:  Supercomputing Center, CNIC, CAS, China
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
	int nPart = 1048576;
	FILE *fp1, *fp2, *fp3;
	char file1[256]="ppn.acc";
	char file2[256]="pp.acc";
	char file3[256]="pp.ori";
	float invg = 1.0;

	if (argc > 1)
	{
		nPart = atoi(argv[1]);
	}
	if (argc > 2)
	{
		strncpy(file2, argv[2], 255);
	}
	if (argc > 3)
	{
		strncpy(file1, argv[3], 255);
	}
	if (argc > 4)
	{
		strncpy(file3, argv[4], 255);
	}
	if (argc > 5)
	{
		invg = 1.0 / atof(argv[5]);
	}


	fp1 = fopen(file1, "r");
	fp2 = fopen(file2, "r");
	fp3 = fopen(file3, "r");

	int n;
	int err[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	for (n=0; n<nPart; n++)
	{
		float x1, y1, z1;
		float x2, y2, z2;
		float dx, dy, dz;
		float x, y, z;
		int m1, m2;
		fscanf(fp1, "%d %f %f %f", &m1, &x1, &y1, &z1);
		fscanf(fp2, "%d %f %f %f", &m2, &x2, &y2, &z2);
		fscanf(fp3, "%d %f %f %f", &m2, &x, &y, &z);
		if (m1 == m2)
		{
			dx = (x2-x1)*invg;
			dy = (y2-y1)*invg;
			dz = (z2-z1)*invg;
			dx = dx>=0 ? dx : -dx;
			dy = dy>=0 ? dy : -dy;
			dz = dz>=0 ? dz : -dz;
			if(dx<=0.0000001f && dy<=0.0000001f && dz<=0.0000001f) err[0]++;
			else if(dx<=0.000001f && dy<=0.000001f && dz<=0.000001f) err[1]++;
			else if(dx<=0.00001f && dy<=0.00001f && dz<=0.00001f) err[2]++;
			else if(dx<=0.0001f && dy<=0.0001f && dz<=0.0001f) err[3]++;
			else
			{
				if(dx<=0.001f && dy<=0.001f && dz<=0.001f)
				{
					err[4]++;
					printf("%d -3  %f/%f  %f/%f  %f/%f  %f  %f  %f\n", m1, x1, x2, y1, y2, z1, z2, x, y, z);
				}
				else if(dx<=0.01f && dy<=0.01f && dz<=0.01f)
				{
					err[5]++;
					printf("%d -2  %f/%f  %f/%f  %f/%f  %f  %f  %f\n", m1, x1, x2, y1, y2, z1, z2, x, y, z);
				}
				else if(dx<=0.1f && dy<=0.1f && dz<=0.1f)
				{
					err[6]++;
					printf("%d -1  %f/%f  %f/%f  %f/%f  %f  %f  %f\n", m1, x1, x2, y1, y2, z1, z2, x, y, z);
				}
				else
				{
					err[7]++;
					printf("%d >-1  %f/%f  %f/%f  %f/%f  %f  %f  %f\n", m1, x1, x2, y1, y2, z1, z2, x, y, z);
				}
			}
		}
	}
	printf ("\n[>-1]%d  [-1]%d  [-2]%d  [-3]%d  [-4]%d  [-5]%d  [-6]%d  [-7]%d\n", err[7], err[6], err[5], err[4], err[3], err[2], err[1], err[0]);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
}
