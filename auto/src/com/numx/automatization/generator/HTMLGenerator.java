package com.numx.automatization.generator;

import java.io.File;
import java.io.FileWriter;
import java.util.Iterator;
import java.util.List;

import org.apache.velocity.VelocityContext;
import org.apache.velocity.context.Context;

import com.numx.automatization.model.Product;
import com.numx.automatization.xml.XMLParser;

public class HTMLGenerator extends VelocityHelper {

	public final static String TEMPLATE_DIR = "com/numx/automatization/generator/doc";
	
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
	}

	public HTMLGenerator() throws Exception{
		super();
		VMs.add(TEMPLATE_DIR);
        initVelocityEngine();
	}

	public void generate(List products, String outDir) throws Exception {
        
        Context ctx = new VelocityContext();
        ctx.put("title", "Product description");
        
    	for (Iterator itProds = products.iterator(); itProds.hasNext();) {
			Product prod = (Product) itProds.next();
	        ctx.put("prod", prod);
	        
	        String outHtmlDir = outDir + "/doc/manual";
	        
	        File f = new File(outHtmlDir, prod.get_productShortNameLower() + "-description.html");
			f.getParentFile().mkdirs();
	        FileWriter fw = new FileWriter(f);
	        velocityEngine.getTemplate(TEMPLATE_DIR +"/signature.vm").merge(ctx, fw);
	        fw.flush();
	        fw.close();
    	}
	}
}
