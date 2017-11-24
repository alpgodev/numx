package com.numx.automatization.generator;

import java.util.List;

import com.numx.automatization.xml.XMLParser;

public class FileGenerator {

	/**
	 * @param args
	 * args[0] is the path of descriptor files
	 * args[1] is the path of root of numX
	 */
	public static void main(String[] args) throws Exception 
	{
		//reading descripor files and getting a list of elements of type Product
		XMLParser parser = new XMLParser();
		
		List products = parser.getProducts(args[0]);
		
		//generating files
		String outDir = args[1];
		System.out.println("Generating files for Functional manual ...");
		new HTMLGenerator().generate(products, outDir);
		
		// API C
		System.out.println("Generating C/C++ API ...");
		new C_API_Generator().generate(products, outDir);
		
		System.out.println("Generating C/C++ solution for Visual Studio ...");
		new C_API_VS2003_Generator().generate(products, outDir);
		
		System.out.println("Generating src for C/C++ tests ...");
		new C_TEST_Generator().generate(products, outDir);
		
		System.out.println("Generating C/C++ test solution for Visual Studio ...");
		new C_TEST_VS2003_Generator().generate(products, outDir);
		
		// API VB
		//System.out.println("Generating Visual Basic API ...");
		//new VB_API_Generator().generate(products, outDir);
		
		//System.out.println("Generating Visual Basic solution for Visual Studio ...");
		//new VB_API_VS2003_Generator().generate(products, outDir);
		
		//System.out.println("Generating src for Visual Basic tests ...");
		//new VB_TEST_Generator().generate(products, outDir);
		
		//System.out.println("Generating Visual Basic test solution for Visual Studio ...");
		//new VB_TEST_VS2003_Generator().generate(products, outDir);
		
		// API Java
		System.out.println("Generating Java API ...");
		new Java_API_Generator().generate(products, outDir);
		
		System.out.println("Generating Java solution for Visual Studio ...");
		new Java_API_VS2003_Generator().generate(products, outDir);
		
		System.out.println("Generating Java test src and build.xml ...");
		new Java_TEST_Generator().generate(products, outDir);
		
		// Scilab API
		//System.out.println("Generating files for Scilab API ...");
		//new ScilabFileGenerator().generate(products, outDir);
		
		// Matlab API
		//System.out.println("Generating files for Matlab API ...");
		//new MatlabFileGenerator().generate(products, outDir);
	}

}
