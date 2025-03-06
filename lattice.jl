
#very stupid way to generate square lattice, but still works
function gen_square_lattice(L)
	function rc2ind(x,y) #row and column to index
		return Int((x-1)*L + y)
	end
	function PBC(x)
		if x>L
			return x-L
		elseif x<1
			return L+x
		else
			return x
		end
	end
	Nc=L*L
	nn=zeros(Int,Nc,4)
	center_pos=zeros(Nc,2)
	#write this function	
	for x in 1:L, y in 1:L
		nn[rc2ind(x,y),:]=[rc2ind(PBC(x+1),y),rc2ind(PBC(x-1),y),rc2ind(x,PBC(y+1)),rc2ind(x,PBC(y-1))]
		center_pos[rc2ind(x,y),1]=x-1
		center_pos[rc2ind(x,y),2]=y-1
	end
	return nn, center_pos
end

#to generate an hexagonal lattice we need to generate a triangular lattice
#to generate a triangular lattice we simply generate a square lattice through the previous function and add diagonal edges
function gen_hex_lattice(L)
	function rc2ind(x,y) #row and column to index
		return (x-1)*L + y
	end
	function PBC(x)
		if x>L
			return x-L
		elseif x<1
			return L+x
		else
			return x
		end
	end
	Nc=L*L
	nn=zeros(Int,Nc,6)
	center_pos=zeros(Nc,2)
	for x in 1:L, y in 1:L
		nn[rc2ind(x,y),:]=[rc2ind(PBC(x+1),y),rc2ind(PBC(x-1),y),rc2ind(x,PBC(y+1)),rc2ind(x,PBC(y-1)),rc2ind(PBC(x-1),PBC(y+1)),rc2ind(PBC(x+1),PBC(y-1))]
		center_pos[rc2ind(x,y),1]=(x-1)+(y%2)*0.5
		center_pos[rc2ind(x,y),2]=(y-1)*sqrt(3)*0.5
		#to have rhombus instead of square use following
		#center_pos[rc2ind(x,y),1]=(y-1)*0.5+(x-1)
		#center_pos[rc2ind(x,y),2]=(y-1)*sqrt(3)*0.5
	end
	return nn,center_pos
end



