#pragma once

class Utilities
{
	public :
		static void setShowDetails(bool);
		static void setShowModules(bool);
		static bool equals(double, double);
		static bool equals(double, double, double);
		static bool equals(int, double[], double[]);
		static bool equals(int, double[], double[], double);
		static bool equals(int, int[], int[]);
		static void printDetailResult(bool, string);	
		static void printModuleResult(bool, string);
		static void printResult(bool, string);
};