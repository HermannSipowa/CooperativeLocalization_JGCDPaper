function    assign( varargin )
%    assign( variable_name, value )  
%    for k = 1:11, assign( sprintf( 'var%d', k ), k ), end
   switch nargin
        case { 2 }
            if  isvarname( varargin{ 1 } )
                Name = varargin{ 1 };
            else
                error( ['poi: First input argument, ', ...
                         inputname(1), ' must be a legal name'] ),
            end
            Value = varargin{ 2 };
        otherwise
            error( 'poi: Wrong number of input arguments' ),
    end
    assignin( 'caller', Name, Value );
end