classdef MyClass < handle
   properties
      a = 0;
      b = [];
   end
   methods
      function obj = myFunc( obj,a, b )
         % Do whatever fmincon wants and set myScalar for it
         obj.a = obj.a + 7;
         obj.b = [ obj.b rand ];
      end
   end
end

% Then in your current code create an object and the function handle, to the function on your object
% myObj = MyClass( );
% f = @myObj.myFunc