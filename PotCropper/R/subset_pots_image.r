
distancetwopoints<-function(point1,point2){
	result<-sqrt((point2[1]-point1[1])^2+(point2[2]-point1[2])^2)
	return(result)
}

readbarcodes_coordinates<-function(dir,file_name){
file_positions<-read.table(paste(dir,file_name,sep=""),header=T)
list_vars<-list()
for (i in 1:nrow(file_positions)){
list_vars[[length(list_vars)+1]]<-list(name=paste("v",i,sep=""),values=c(file_positions$row[i],file_positions$column[i]))
}
return(list_vars)
}

combn_distances_barcodes<-function(list_vars){
	combi_individuals<-combn(sapply(list_vars,"[[","name"),2)
	distan_comb<-list()
	for(i in 1:ncol(combi_individuals)){
		part_1<-list_vars[which(sapply(list_vars,"[[","name")==combi_individuals[,i][1])]
		part_2<-list_vars[which(sapply(list_vars,"[[","name")==combi_individuals[,i][2])]
		p1<-as.vector(sapply(part_1,"[[","values"))
		p2<-as.vector(sapply(part_2,"[[","values"))
		length<-distancetwopoints(p1,p2)
		distan_comb[[length(distan_comb)+1]]<-list(name=paste(combi_individuals[,i],collapse="_"),values=length)
	}
	return(distan_comb)
}

order_list_by_numeric<-function(list){
	distances <- sapply(list,"[[","values")
	identifiers <- sapply(list,"[[","name")
	atomic_elements<-strsplit(identifiers[order(distances,decreasing=T)],"_")
	return(atomic_elements)
}
determine_h1_h2_h3<-function(dir,file){
indiv_vars<-readbarcodes_coordinates(dir,file)
indiv_variables<-sapply(indiv_vars,"[[","name")
atomic_elements<-order_list_by_numeric(combn_distances_barcodes(indiv_vars))
h1<-intersect(atomic_elements[[1]],atomic_elements[[2]])
h2<-atomic_elements[[1]][atomic_elements[[1]]!=h1]
h3<-indiv_variables[indiv_variables!=h1 & indiv_variables !=h2]
h1<-list(sapply(indiv_vars[which(sapply(indiv_vars,"[[","name")==h1)],"[[","values"))
h2<-list(sapply(indiv_vars[which(sapply(indiv_vars,"[[","name")==h2)],"[[","values"))
h3<-list(sapply(indiv_vars[which(sapply(indiv_vars,"[[","name")==h3)],"[[","values"))
final<-c(h1,h2,h3)
return(final)
}

slopetwopoints<-function(point1,point2){
	result<-(point2[1]-point1[1])/(point2[2]-point1[2])
	return(result)
}

determinationofb<-function(point,slope){
	result<-point[1]-(slope*point[2])
	return(result)
}
lineform<-function(point1,point2){
	slope<-slopetwopoints(point1,point2)
	b<-determinationofb(point1,slope)
	return(polynomial(c(b,slope)))

}
evalpolynom<-function(polynom,value){
	toevaluate<-as.function(polynom)
	result<-toevaluate(value)
	return(result)
}

substitute_variable_polynom<-function(poly_original,poly_insert){
	final_polynom<-poly_original[1]
	for (i in 2:length(poly_original)){
		final_polynom<-(final_polynom)+((poly_original[i])*(as.polynomial(poly_insert))^(i-1))
	}
	

	return(final_polynom)
}

quadraticform<-function(variable_coefficient,zero_terminus){
first<-polynomial(c(zero_terminus,variable_coefficient))
second<-first^2
return(second)
}

find_xs<-function(point,polynom,distance){
	sols_x<-solve(quadraticform(1,-point[2])+substitute_variable_polynom(quadraticform(1,-point[1]),polynom)- (distance)^2)
	return(sols_x)
}

find_ys<-function(sols_x,line_pol){
	sols<-list()
	for(i in 1:length(sols_x)){
		sols[[i]]<-c(evalpolynom(line_pol,sols_x[i]),sols_x[i])
	}
	return(sols)
}


getsols<-function(point,polynom,distance){
return(find_ys(find_xs(point,polynom,distance),polynom))
}

returncloserpoint<-function(sols_xy,ref){
		a<-distancetwopoints(unlist(sols_xy[[1]]),ref)
		b<-distancetwopoints(unlist(sols_xy[[2]]),ref)
		if(a<b){
			return(1)
		}
		else{
			return(2)
		}	
}


findmidpoint<-function(h1,h2){
	polynom<-lineform(h1,h2)
	sols_xy<-getsols(h1,polynom,(distancetwopoints(h1,h2)/2))
	result<-sols_xy[returncloserpoint(sols_xy,h2)]
	return(result)
}

straightline<-function(point,slope){
	b<-determinationofb(point,slope)
	return(polynomial(c(b,slope)))
}

intersect_twolines<-function(pol1,pol2){
	x<-solve(polynomial(c(pol1[1]+(pol2[1]*-1),pol1[2]+(pol2[2]*-1))))
	y<-evalpolynom(pol1,x)
	return(c(y,x))
}

correct_h1_h2_h3<-function(h1,h2,h3,position){
	final<-list()
	
	sols_h1<-getsols(h1,straightline(h1,(1/slopetwopoints(h1,h3))*-1),130)
	sols_h3<-getsols(h3,straightline(h3,(1/slopetwopoints(h1,h3))*-1),130)
	sols_h2<-getsols(h2,straightline(h2,(1/slopetwopoints(h1,h3))*-1),130)
	if (position == "outside"){
		final_h1<-sols_h1[returncloserpoint(sols_h1,h2)]
		final_h3<-sols_h3[returncloserpoint(sols_h3,h2)]
		final_h2<-sols_h2[returncloserpoint(sols_h2,h3)]

	}
	else{
		final_h1<-sols_h1[-returncloserpoint(sols_h1,h2)]
		final_h3<-sols_h3[-returncloserpoint(sols_h3,h2)]
		final_h2<-sols_h2[-returncloserpoint(sols_h2,h3)]
	}
	final[1]<-list(final_h1)
	final[2]<-list(final_h2)
	final[3]<-list(final_h3)
	return(final)
}



findc1andc2<-function(h1,h2,h3){
	final<-list()
	c2<-intersect_twolines(straightline(h1,slopetwopoints(h1,h3)),straightline(h2,(1/slopetwopoints(h1,h3))*-1))
	c1<-intersect_twolines(straightline(h1,(1/slopetwopoints(h1,h3)*-1)),straightline(h2,slopetwopoints(h1,h3)))
	final[2]<-list(c2)
	final[1]<-list(c1)
	return(final)

}


matrix_points<-function(h1,c1,c2,midpoint){
slope_h1_c2<-slopetwopoints(h1,c2)
p_h1_c2<-straightline(h1,slope_h1_c2)
dis_h1_c2<-distancetwopoints(h1,c2)
prop_dis_h1_c2<-(1/3)*dis_h1_c2
slope_h1_c1<-slopetwopoints(h1,c1)
p_h1_c1<-straightline(h1,slope_h1_c1)
dis_h1_c1<-distancetwopoints(h1,c1)
prop_dis_h1_c1<-(1/3)*dis_h1_c1
fin<-matrix(data=NA,16,4,byrow=T)
rnames<-c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14","p15","p16")
cnames<-c("Row_Coordinate","Column_Coordinate","Orientarion_Column","Orientation_Row")
k=1
for (i in 0:3){

	if((dis_h1_c2 - prop_dis_h1_c2*i) >=0){
		if (i ==0){
			refpoint<-h1
			fin[k,1]<-refpoint[1]
			fin[k,2]<-refpoint[2]
			fin[k,3]<-p_h1_c2[2]
			fin[k,4]<-p_h1_c1[2]
			k=k+1
		}
		else{
		sols<-getsols(h1,p_h1_c2,(prop_dis_h1_c2*i)+0)
		points_relation_square<-returncloserpoint(sols,unlist(midpoint))
		refpoint<-sols[[points_relation_square]]
		fin[k,1]<-refpoint[1]
		fin[k,2]<-refpoint[2]
		fin[k,3]<-p_h1_c2[2]
		fin[k,4]<-p_h1_c1[2]
		k=k+1
		}
		for (j in 1:3){
		if((dis_h1_c1 - prop_dis_h1_c1*j) >=0){
		sols_2<-getsols(refpoint,straightline(refpoint,slopetwopoints(h1,c1)),(prop_dis_h1_c1*j)+0)
		#sols_coordinates_2<-sols_into_coordinates(sols_2)
		points_relation_square_2<-returncloserpoint(sols_2,unlist(midpoint))
		#points_relation_square_2<-point.in.polygon(sols_coordinates_2[[1]],sols_coordinates_2[[2]],square_x,square_y)
		#Esta parte muere en las orillas porque considero la liena entre c1 y h2, debo poner
		#print(sols_2[[points_relation_square_2]])
		refpoint_2<-sols_2[[points_relation_square_2]]
		fin[k,1]<-refpoint_2[1]
		fin[k,2]<-refpoint_2[2]
		fin[k,3]<-p_h1_c2[2]
		fin[k,4]<-p_h1_c1[2]
		k=k+1
			}
		}
	}
}
rownames(fin)<-rnames
colnames(fin)<-cnames
return(fin)
}


pots_from_point_matrix<-function(image,matrix,i,size=length){
	point<-c(matrix[i,1],matrix[i,2])
	slope_y<-matrix[i,4]
	slope_x<-matrix[i,3]

	sols_row<-getsols(point,straightline(point,slope_y),size)
	#points_relation_square<-sols_row[[returncloserpoint(sols_row,point)]]
	sols_col<-getsols(point,straightline(point,slope_x),size)
	#points_relation_square<-sols_col[[returncloserpoint(sols_col,point)]]
	if(abs(sols_row[[2]][1]-sols_row[[1]][1])<size){
	section<-image[sols_row[[1]][2]:sols_row[[2]][2],sols_col[[1]][1]:sols_col[[2]][1],]

	}
	else{
	section<-image[sols_col[[1]][2]:sols_col[[2]][2],sols_row[[1]][1]:sols_row[[2]][1],]
	}
	return(section)

}



pmatrix_from_image<-function(dir,file,position){
	hs<-determine_h1_h2_h3(dir,file)
	temp_h1<-unlist(hs[1])
	temp_h2<-unlist(hs[2])
	temp_h3<-unlist(hs[3])
	hs_corrected<-correct_h1_h2_h3(temp_h1,temp_h2,temp_h3,position)
	h1<-unlist(hs_corrected[1])
	h2<-unlist(hs_corrected[2])
	h3<-unlist(hs_corrected[3])
	midpoint<-findmidpoint(h1,h2)
	cs<-findc1andc2(h1,h2,h3)
	c1<-unlist(cs[1])
	c2<-unlist(cs[2]) 	
	mpoints<-matrix_points(h1,c1,c2,midpoint)
	return(mpoints)
}

#' Subset pots v_0.2
#'
#' This functions takes a single image or all the images contained in a directory and subset all the pots contained in the image
#' By default this function only takes JPG images (case insensitive).
#' It is mandatory that the last character of  dir and outputdir is "/" as in "/path1/path2/path3/"
#' The default dir and outputdir is "./"
#' User should specify the position of the barcode modyfing the parameter -barcode_position- Default outside
#' -pots_positions- defined which specific pots of the image the user want to subset. Default 1-16
#' -size- defined the length of the square surrounding the pot to subser. Default 350
#' Examples of Usage
#'
#' Single image in ./ : subset_pots("my.image.jpg")
#' Single image in other path :subset_pots("/path/file.jpg")
#' Full directory :subset_pots(dir="path/")
#' Full directory ./ : subset_pots()
subset_pots<-function(file_image=NA,dir="./",outdir="./",file_txt=paste(unlist(strsplit(paste(file_image),"[.]"))[1],"txt",sep="."),barcode_position="outside",pots_positions=c(1:16),size=325){
	if (is.na(file_image)){
		imfiles<-list.files(path=dir)[c(grep(".JPG",list.files(path = dir),ignore.case=T))]
		for (i in 1:length(imfiles)){
			imagen<-readImage(paste(dir,imfiles[i],sep=""))
			file_id<-unlist(strsplit(paste(imfiles[i]),"[.]"))[1]
			file_txt<-paste(unlist(strsplit(paste(imfiles[i]),"[.]"))[1],"txt",sep=".")
			for(i in 1:length(pots_positions)){
				section<-pots_from_point_matrix(imagen,pmatrix_from_image(dir,file_txt,barcode_position),i,size)
				writeImage(section,paste(paste(outdir,paste(file_id,"pot",i,sep="_"),sep=""),"jpeg",sep="."))
			}
		}
	}
	else if (!is.na(file_image)){
		parts<-unlist(strsplit(file_image,"/"))
		if(length(parts) > 1){
			rfile<-tail(parts,n=1)
			ldir<-parts[parts!=tail(parts,n=1)]
			fdir<-paste(ldir,collapse="/")
			fdir<-paste(fdir,"/",sep="")
			imagen<-readImage(paste(fdir,rfile,sep="/"))
			file_id<-unlist(strsplit(paste(rfile),"[.]"))[1]
			file_txt<-paste(unlist(strsplit(tail(unlist(strsplit(file_image,"/")),n=1),"[.]"))[1],"txt",sep=".")
			for(i in 1:length(pots_positions)){
				section<-pots_from_point_matrix(imagen,pmatrix_from_image(fdir,file_txt,barcode_position),i,size)
				writeImage(section,paste(paste(outdir,paste(file_id,"pot",i,sep="_"),sep=""),"jpeg",sep="."))

		}
	}
		else{
		imagen<-readImage(paste(dir,file_image,sep=""))
		file_id<-unlist(strsplit(paste(file_image),"[.]"))[1]
			for(i in 1:length(pots_positions)){
				section<-pots_from_point_matrix(imagen,pmatrix_from_image(dir,file_txt,barcode_position),i,size)
				writeImage(section,paste(paste(outdir,paste(file_id,"pot",i,sep="_"),sep=""),"jpeg",sep="."))
	
	}
	
}
}
}
