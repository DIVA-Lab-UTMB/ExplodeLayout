#' Explode the old coordinates using cluster membership
#'
#' Takes the nodelist of a network and return an updated nodelist with the exploded coordinates.
#'
#' @param nodelist A nodelist is a dataframe including at least 4 columns: Label, X, Y, Cluster. Label is the node ID. They must be unique character strings. X and Y are the original coordinates. Cluster is the node clustering membership. They must be integers ranging from 1 to the number of clusters.
#' @param radius The radius of the used to explode the clusters, which defaults to 1.
#'
#' @return A nodelist with new node layout coordinates generated from the algorithm. This is a new dataframe with 4 colums: Label, newX, newY, Cluster.
#'
#' @examples
#' exploded_coords=explode_coordinates(example_nodelist,radius=1.2)
explode_coordinates=function(nodelist,radius=1){
  #nodelist=data.frame(Label,X,Y,Cluster)

  nclust=length(unique(nodelist$Cluster)) # number of clusters
  net_center=list('X'=mean(nodelist$X),'Y'=mean(nodelist$Y))
  net_radius=sqrt(sum((nodelist$X-net_center$X)^2+(nodelist$Y-net_center$Y)^2)/length(nodelist$Label)) # standard deviation of nodes' distance to center
  radius_unit=nclust*net_radius/pi # heuristic radius

  clu_centers=as.data.frame(matrix(0,nrow=nclust,ncol=3))
  colnames(clu_centers)=c('X','Y','angle')
  for(i in seq(nclust)){
    clu_nl=nodelist[nodelist$Cluster==i,]
    clu_center=list('X'=mean(clu_nl$X),'Y'=mean(clu_nl$Y))
    clu_centers[i,'X']=clu_center$X
    clu_centers[i,'Y']=clu_center$Y
    clu_centers[i,'angle']=atan2(clu_center$Y-net_center$Y,clu_center$X-net_center$X) # original angle -pi to pi
  }
  clu_centers$rank=rank(clu_centers$angle,ties.method = 'first')
  clu_centers$new_angle=(clu_centers$rank-1)/nclust*2*pi-pi # uniform distribution from -pi to pi
  clu_centers$angle_diff=clu_centers$new_angle-clu_centers$angle # rotation
  clu_centers$newX=radius_unit*radius*cos(clu_centers$new_angle)+net_center$X # new X for cluster center
  clu_centers$newY=radius_unit*radius*sin(clu_centers$new_angle)+net_center$Y # new Y for cluster center

  nodelist$shiftedX=0
  nodelist$shiftedY=0
  nodelist$rotatedX=0
  nodelist$rotatedY=0
  for(i in seq(nclust)){
    clu_nodes=(nodelist$Cluster==i)
    nodelist[clu_nodes,'shiftedX']=nodelist[clu_nodes,'X']+clu_centers[i,'newX']-clu_centers[i,'X']
    nodelist[clu_nodes,'shiftedY']=nodelist[clu_nodes,'Y']+clu_centers[i,'newY']-clu_centers[i,'Y']
    nodelist[clu_nodes,'rotatedX']=cos(clu_centers[i,'angle_diff'])*(nodelist[clu_nodes,'shiftedX']-clu_centers[i,'newX'])-sin(clu_centers[i,'angle_diff'])*(nodelist[clu_nodes,'shiftedY']-clu_centers[i,'newY'])+clu_centers[i,'newX']
    nodelist[clu_nodes,'rotatedY']=sin(clu_centers[i,'angle_diff'])*(nodelist[clu_nodes,'shiftedX']-clu_centers[i,'newX'])+cos(clu_centers[i,'angle_diff'])*(nodelist[clu_nodes,'shiftedY']-clu_centers[i,'newY'])+clu_centers[i,'newY']
  }
  nodelist$newX=nodelist$rotatedX
  nodelist$newY=nodelist$rotatedY
  return(nodelist[c('Label','newX','newY','Cluster')])
}

#' Generate edgelist from incidence matrix of a bipartite network
#'
#' @param incidence_matrix A matrix where row names and column names are the node ID of a bipartite network. An element of the i-th row and j-th column of the matrix is 0 if node on row i is not connected to node on column j, and edge weight if they are connected.
#'
#' @return A dataframe with 3 columns: nodesR, nodesC, values.
#'
#' @examples
#' example_edgelist=get_edgelist_from_incidmat(example_incidmat)
get_edgelist_from_incidmat=function(incidence_matrix){
  el=utils::stack(incidence_matrix) # a dataframe of two columns: values, ind
  el$nodesR=rownames(incidence_matrix)
  el$nodesC=as.character(el$ind)
  return(el[el$values!=0,c('nodesR','nodesC','values')])
}

#' Format a nodelist for plotting using the new coordinates found by ExplodeLayout
#'
#' @param nodelist A dataframe including at least 3 columns: Label, Cluster, Entity. Entity is either 1 or 2, indicating which part of the bipartite network a node is in.
#' @param new_coordinates A dataframe including at least 3 columns: Label, newX, newY.
#'
#' @return A dataframe including 5 columns: Label, X, Y, Color, baseShape.
#'
#' @examples
#' exploded_coords=explode_coordinates(example_nodelist,radius=1.2)
#' plotting_nodelist=get_nodelist_for_plotting(example_nodelist,exploded_coords)
get_nodelist_for_plotting=function(nodelist,new_coordinates){
  nl=nodelist
  rownames(nl)=nl$Label
  rownames(new_coordinates)=new_coordinates$Label
  new_coordinates=new_coordinates[rownames(nl),]
  nl$X=new_coordinates$newX
  nl$Y=new_coordinates$newY
  #colorList=c('red','orange','green','blue')
  colorList=c("#C9197A", "#94EA18", "#0AF3EE", "#808080", "#E5F115", "#8825F9", "#CE9A89",  "#fd991c", "black", "#6b8e23", "#C0C4E0", "#669999", "blue", "pink") # old explode layout
  nl$Color=colorList[nl$Cluster]
  nl$baseShape=nl$Entity
  nl[nl$Entity==1,'baseShape']=21
  nl[nl$Entity==2,'baseShape']=24
  nl$igraphShape=nl$Entity
  nl[nl$Entity==1,'igraphShape']='circle'
  nl[nl$Entity==2,'igraphShape']='triangle'
  nl$PrintLabel=nl$Label
  nl[nl$Entity==1,'PrintLabel']=NA
  return(nl[,c('Label','X','Y','Color','baseShape')])
}

#' Explode the old coordinates using cluster membership and generate nodelist for plotting
#'
#' Takes the nodelist of a network and return an updated nodelist with the exploded coordinates based on radius, node color based on cluster, and node shape based on entity.
#'
#' @param nodelist A nodelist is a dataframe including at least 5 columns: Label, X, Y, Cluster, Entity. Label is the node ID. They must be unique character strings. X and Y are the original coordinates. Cluster is the node clustering membership. They must be integers ranging from 1 to the number of clusters. Entity indicates which part of the bipartite network a node belongs to. (Can be either 1 or 2.)
#' @param radius The explode radius of the projecting circle. Default to 1.
#'
#' @return A new nodelist with exploded coordinates, which is a dataframe including 5 columns: Label, X, Y, Color, baseShape.
#' @export
#'
#' @examples
#' exploded_nodelist=get_explode_nodelist(example_nodelist,radius=1.2)
get_explode_nodelist=function(nodelist,radius=1){
  exploded_coords=explode_coordinates(nodelist,radius=1.2)
  plotting_nodelist=get_nodelist_for_plotting(nodelist,exploded_coords)
  return(plotting_nodelist)
}


#' Plot bipartite network given node list (label, coordinates, shape, color) and incidence matrix.
#'
#' @param nodelist A dataframe of at least 5 columns: Label, X, Y, Color, baseShape.
#' @param incidence_matrix A matrix where row names and column names are the node ID of a bipartite network. An element of the i-th row and j-th column of the matrix is 0 if node on row i is not connected to node on column j, and edge weight if they are connected.
#'
#' @return a ggplot2 object p which can be shown using print(p).
#' @export
#'
#' @examples
#' exploded_nodelist=get_explode_nodelist(example_nodelist,radius=1.2)
#' p=plot_binet_ggplot2(exploded_nodelist,example_incidmat)
#' print(p)
plot_binet_ggplot2=function(nodelist,incidence_matrix){
  rownames(nodelist)=nodelist$Label
  el=get_edgelist_from_incidmat(incidence_matrix)
  el$x0=nodelist[el$nodesR,'X']
  el$y0=nodelist[el$nodesR,'Y']
  el$x1=nodelist[el$nodesC,'X']
  el$y1=nodelist[el$nodesC,'Y']

  p=ggplot2::ggplot()
  #p=p + ggplot2::geom_segment(ggplot2::aes(x = el$x0, y = el$y0, xend = el$x1, yend = el$y1), data = el, size = 0.01, colour = "#888888")
  p=p + ggplot2::geom_segment(ggplot2::aes(x = el$x0, y = el$y0, xend = el$x1, yend = el$y1), size = 0.01, colour = "#888888")
  #p=p + ggplot2::geom_point(ggplot2::aes(nodelist$X, nodelist$Y),colour = "black",  fill = nodelist$Color, shape = nodelist$baseShape, size = 3, data = nodelist)
  p=p + ggplot2::geom_point(ggplot2::aes(nodelist$X, nodelist$Y),colour = "black",  fill = nodelist$Color, shape = nodelist$baseShape, size = 3)
  p=p + ggplot2::theme(panel.background = ggplot2::element_blank()) + ggplot2::theme(legend.position = "none")
  p=p + ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank())
  p=p + ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank())
  return(p)
}
