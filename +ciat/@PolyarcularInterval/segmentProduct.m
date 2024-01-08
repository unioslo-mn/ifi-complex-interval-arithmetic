function outObj = segmentProduct(obj1, obj2)


	% Check input class
    mustBeA(obj1,["ciat.Edge","ciat.Arc"]);
    mustBeA(obj2,["ciat.Edge","ciat.Arc"]);
    


	if type(obj1) = ciat.Edge
		if type(obj1) = ciat.Edge

		else

		end
	else
		if type(obj1) = ciat.Edge

		else

		end
	end

end

function out = edgeTimesEdge(edge1,edge2)

	if edge1.CurveParam == 0 || edge2.CurveParam == 0

	else

	end

end

function out = edgeTimesArc(edge,arc)

	if edge.CurveParam == 0

	elseif arc.CurveParam == 0

	else

	end

end

function out = arcTimesArc(arc1,arc)


	if arc1.CurveParam == 0 || arc2.CurveParam == 0

	else

	end

end
