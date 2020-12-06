#' Load STL file
#' 
#' This is a wrapper over RGL's readSTL function
#' 
#' @param con A connection or filename
#' @return A matrix with x, y, and z points
#' @export
stl_load <- function(con) {
	stl <- rgl::readSTL(con=con, plot=FALSE)
	colnames(stl) <- c("x", "y", "z")
	return(stl)
}

#' Get dimensions of the object
#' 
#' @param stl object from stl_load()
#' @return matrix with min and max values for x, y, and z
#' @export
stl_dimensions <- function(stl) {
	return(apply(stl, 2, range)) 
}

# #' Make rectangular array
# #' 
# #' The object may only be defined for certain points: if it is a cone, it may have a circular base, but we want to cut a rectangular space around it. This will pad that space (x and y only). You can choose what units to use: "raw" will add in the dimensions of the object (if it is going from -40 to 40 mm, then padding=1, padding_units="raw" will mean add 1 mm on each side). "percent" will increase it by a percentage of the relevant axis: padding=1, padding_units="percent" will make it 1% wider on each side.
# #' 
# #' The zmax affects whether the new points are placed at the top of the object (zmax), as you might want for the Grand Canyon, or at the minimum, as you might want if doing a mountain.
# #' 
# #' @param stl object from stl_load()
# #' @param padding how much to pad it by
# #' @param padding_units 'raw' units (the units of the object) or 'percent'
# #' @param zmax If TRUE, the new points are at the maximum value of Z; if FALSE, the minimum.
# #' @return A matrix with x, y, and z points
# stl_pad <- function(stl, padding=0, padding_units='raw', zmax=TRUE) {

# }

#' Regularize the object
#' 
#' The STL object might only be defined over a certain area: a model of a mountain might stop at the irregular base. The points might also not be on a regular grid, which makes the analysis in R difficult. So we will use a technique called kriging to smoothly interpolate the points. The question is how fine to make the final grid: does it have as many points as the original points, or is it finer? By default it will have 5 points between every given point on each dimension. If it's too large, consider using meshlab to scale it down well before loading here.
#' 
#' zero_position resets the x and y positions so that (0,0) can be the bottom left of the structure ("bottomleft"), the center ("center"), or kept as is ("current")
#'
#' STL files might need to be rescaled. If  set_max_dimension is set to a number, the x, y, and z dimensions are set so that the largest width or length is this value and the other dimensions are scaled correctly.
#'
#' @param stl object from stl_load()
#' @param fineness how many points to make for each point per dimension in the original stl
#' @param nmax how many nearby points to use for kriging. Inf to use all points.
#' @param zero_position_xy can be "bottomleft", "center", or "current".
#' @param zero_position_z can be "top", "bottom", or "current".
#' @param set_max_dimension if not NULL, rescales the object so the maximum width or length is this size
#' @param verbose if TRUE, report on progress
#' @return A matrix with x, y, and z points; x and y are evenly spaced grid
#' @export
#' @examples
#' cat_paw <- stl_load("~/Downloads/catsnowscaled.stl")
#' cat_paw_reg <- stl_regularize(cat_paw)
stl_regularize <- function(stl, fineness=10, nmax=20, zero_position_xy="bottomleft", zero_position_z="top", max_dimension=NULL, verbose=TRUE) {


	if(verbose) {
		print(paste0("The stl is ", xwidth, " units wide (left to right, x), ", ywidth, " units tall (forward and back, y), and ", zwidth, " units deep (up and down into the material, z). Typically these units are mm, but it may vary with CNC."))
	}
	if(zero_position_xy=="bottomleft") {
		stl[,"x"] <- stl[,"x"]-min(stl[,"x"])
		stl[,"y"] <- stl[,"y"]-min(stl[,"y"])
	}
	if(zero_position_xy=="center") {
		stl[,"x"] <- 0.5*xwidth + stl[,"x"]-min(stl[,"x"])
		stl[,"y"] <- 0.5*ywidth + stl[,"y"]-min(stl[,"y"])
	}
	if(zero_position_z=="top") {
		stl[,"z"] <- stl[,"z"]-max(stl[,"z"])
	}
	if(zero_position_z=="bottom") {
		stl[,"z"] <- stl[,"z"]-min(stl[,"z"])
	}


	ranges <- as.data.frame(stl_dimensions(stl))
	xwidth <- abs(diff(ranges$x))
	ywidth <- abs(diff(ranges$y))
	zwidth <- abs(diff(ranges$z))
	maxwidth <- max(xwidth, ywidth) # so we do a grid with same resolution for x and y

	if(!is.null(set_max_dimension)) {
		scaling_factor <- set_max_dimension/maxwidth
		stl[,"x"] <- stl[,"x"] * scaling_factor
		stl[,"y"] <- stl[,"y"] * scaling_factor
		stl[,"z"] <- stl[,"z"] * scaling_factor
		ranges <- as.data.frame(stl_dimensions(stl))
		xwidth <- abs(diff(ranges$x))
		ywidth <- abs(diff(ranges$y))
		zwidth <- abs(diff(ranges$z))
		maxwidth <- max(xwidth, ywidth) # so we do a grid with same resolution for x and y
	}

	minratio <- min(xwidth, ywidth)/maxwidth
	original_npoints_per_dim <- round(sqrt(fineness*nrow(stl)/minratio))
	
	(fineness*max(length(unique(stl[,"x"])), length(unique(stl[,"y"]))))
	# we might hit memory limits
	new_array <- NULL
	npoints_per_dim <- original_npoints_per_dim
	stl_spatial <- sp::SpatialPoints(stl)
	while(class(new_array)!="SpatialPoints") {
		if(verbose) {
			print(paste0("Now trying to create new array with ", round(npoints_per_dim*xwidth/maxwidth), " points per width (x) and ", round(npoints_per_dim*ywidth/maxwidth), " points per height (y) for ", round(npoints_per_dim*xwidth/maxwidth) * round(npoints_per_dim*ywidth/maxwidth), " points total. The original object has ", nrow(stl), " points total, but these are typically not evenly spaced; the new matrix will be ", round(round(npoints_per_dim*xwidth/maxwidth) * round(npoints_per_dim*ywidth/maxwidth) / nrow(stl),1), " times bigger."))
		}
		try(new_array <- sp::SpatialPoints(expand.grid(x=seq(from=min(ranges$x), to=max(ranges$x), length.out=round(npoints_per_dim*xwidth/maxwidth)), y=seq(from=min(ranges$y), to=max(ranges$y), length.out=round(npoints_per_dim*ywidth/maxwidth)))))
		npoints_per_dim <- 0.1*npoints_per_dim
	}
	if(original_npoints_per_dim != 10*npoints_per_dim) {
		warning("The array with the desired fineness was too large for R, so the resolution was decreased")
	}
	print("Now doing interpolation (kriging) to get the interpolated depths")
	new_stl <- as.matrix(as.data.frame(gstat::krige(z~1, stl_spatial, newdata=new_array, nmax=nmax))[,1:3])
	colnames(new_stl) <- c("x", "y", "z")
	return(new_stl)
}


#' Creates the GCODE for a regularized loaded stl object
#' 
#' The defaults work for my machine (a MySweety 3018) but you should check this for yours. You can look at the final code; the site https://www.cnccookbook.com/g-code-m-code-reference-list-cnc-mills/ will explain what codes mean. The important things to check are the three speeds, the units you want to use, and the stepover_height. spin_speed is how quickly to rotate the spindle. horizontal_speed is how quickly to move left and right and forward and back while cutting through material: too fast and it might be too hard for your machine (or snap your bit), too slow and the run takes longer. Vertical is the same but for moving into and out of the material. stepover_height is how much clearance to leave above obstacles -- make sure your spindle can go up at least this much above its starting position. stepdown_depth is how much lower to go into the material each pass.
#' 
#' Remember, this code is not guaranteed to work! And it is code telling an expensive robot how to move very quickly moving bits of metal across resisting surfaces. If it goes wrong, it can go VERY wrong: damaging your machine, damaging the materials, setting things on fire, hurting people. Use at your own risk, and I strongly suggesting seeking out better options (like FreeCAD) if they are available to you. The settings seem to work for my machine, but they could be too fast or slow for yours, make sure you don't try to engrave something too large, etc. You can run your code at https://nraynaud.github.io/webgcode/ to see what it will look like and how long it will take.
#' 
#' The algorithm I use here is very basic -- essentially moving back and forth like a farmer plowing their field. There are far more clever algorithms like spiraling, climbing cuts, etc. This also does not focus on the tool size: basically, it moves the bottom tip of the tool back and forth along each line specified in the object, moving a little bit downward each time. 
#' 
#' @param stl object from stl_load(), ideally after running through stl_regularize
#' @param gcode_file path to use to store gcode output
#' @param spin_speed how fast to rotate the spindle
#' @param horizontal_speed how fast to move the spindle when cutting horizontally
#' @param vertical_speed how fast to move the spindle when moving down or up
#' @param stepover_height how high to move the spindle over an obstacle
#' @param stepdown_depth how much lower to move the spindle each pass
#' @param unit "inches" or "mm"
#' @param unit_precision how many digits to use after the decimal point
#' @param verbose if TRUE, print details
#' @export
#' @examples
#' cat_paw <- stl_load("~/Downloads/catsnowscaled.stl")
#' cat_paw_reg <- stl_regularize(cat_paw)
#' stl_generate_gcode(cat_paw_reg, gcode_file="~/Downloads/cat.nc")
stl_generate_gcode <- function(stl, gcode_file='gcode.nc', spin_speed=12000, horizontal_speed=30, vertical_speed=9, stepover_height = 0.2, stepdown_depth=0.02, unit="inches", unit_precision=5, verbose=TRUE) {
	ranges <- as.data.frame(stl_dimensions(stl))
	stl <- round(stl, unit_precision)
	xwidth <- abs(diff(ranges$x))
	ywidth <- abs(diff(ranges$y))
	zwidth <- abs(diff(ranges$z))
	if(verbose) {
		print(paste0("The stl is ", round(xwidth,3), " ", unit, " wide (left to right, x), ", round(ywidth,3), " ", unit, " tall (forward and back, y), and ", round(zwidth,3), " ", unit, " deep (up and down into the material, z)."))
	}
	cat(ifelse(unit=="mm", "G21\n", "G20\n"), file=gcode_file, append=FALSE) #set the units used
	cat(paste0("M3 S", spin_speed, "\n"), file=gcode_file, append=TRUE) #start spindle turning clockwise at specified speed
	cat("G90\n", file=gcode_file, append=TRUE) #specify absolute coding
	
	depth_passes <- seq(from=0, to=min(stl[,"z"]), by=-1*abs(stepdown_depth))[-1] #how deep each pass should be, but eliminating the first pass at zero depth
	x_positions <- sort(unique(stl[,"x"]), decreasing=FALSE)
	y_positions <- sort(unique(stl[,"y"]), decreasing=FALSE)
	for (depth_index in seq_along(depth_passes)) {
		starting_y <- min(stl[,"y"])
		starting_x <- min(stl[,"x"])
		cat(paste0("G1 Z", stepover_height, " F", vertical_speed, "\n"),  file=gcode_file, append=TRUE) #let's slowly go up to clear the piece
		cat(paste0("G0 X", starting_x, " Y", starting_y, "\n"), file=gcode_file, append=TRUE) #shoot over to the starting position
		desired_z <- depth_passes[depth_index]
		for (x_index in seq_along(x_positions)) {
			x_carve_position <- x_positions[x_index]
			stl_slice <- stl[stl[,"x"]==x_positions[x_index],]
			stl_slice <- stl_slice[order(stl_slice[,"y"], decreasing=FALSE),]
			valid_stl_indices <- (stl_slice[,"z"]<=desired_z)
			if(x_index%%2==0) { #do odd passes in one direction, even passes in the other, to result in less travel
				valid_y_indices <- rev(valid_y_indices)
			}

# now from this vector figure out where the start and stop are: go to the first start, move down, move over to the first stop, pull up, move to the next start, etc. Worry about an odd number of segments (there could be one position that is a start and a stop)
			if(length(valid_y_indices)>0) {
				cat(paste0("G0 X", x_carve_position, " Y", y_positions[valid_y_indices[1]], "\n"), file=gcode_file, append=TRUE) #shoot over to the first position where we will plunge for this row
				cat(paste0("G1 Z", desired_z, " F", vertical_speed, "\n"),  file=gcode_file, append=TRUE) #plunge to start the cut
				valid_y_indices <- valid_y_indices[-1] #we took care of the first point, now repeat
				while(length(valid_y_indices)>0) {
					cat(paste0("G1 X", x_carve_position, " Y", y_positions[valid_y_indices[1]], " F", horizontal_speed, "\n"), file=gcode_file, append=TRUE) #shoot over to the first position where we will plunge for this row
					valid_y_indices <- valid_y_indices[-1] #we now moved over to the last point before having to go up
					cat(paste0("G1 Z", stepover_height, " F", vertical_speed, "\n"),  file=gcode_file, append=TRUE) #let's slowly go up to clear the piece
					if(length(valid_y_indices)>0) {
						cat(paste0("G0 X", x_carve_position, " Y", y_positions[valid_y_indices[1]], "\n"), file=gcode_file, append=TRUE) #shoot over to the first position where we will plunge for this row
						cat(paste0("G1 Z", desired_z, " F", vertical_speed, "\n"),  file=gcode_file, append=TRUE) #plunge to start the cut
						valid_y_indices <- valid_y_indices[-1] #we now handled the next point
					}
				}
			}
			cat(paste0("G1 Z", stepover_height, " F", vertical_speed, "\n"),  file=gcode_file, append=TRUE) #let's get out of the way
		} # end x positions loop for this depth

		
	} # end depth_passes loop
	cat(paste0("G1 Z", stepover_height, " F", vertical_speed, "\n"),  file=gcode_file, append=TRUE) #let's get out of the way, though should already be clear
	cat(paste0("G0 X", starting_x, " Y", starting_y, "\n"), file=gcode_file, append=TRUE) #shoot over to the starting position
	cat(paste0("M05"), file=gcode_file, append=TRUE) #shoot over to the starting position
	print(paste0("Done, saved to ", gcode_file))
}