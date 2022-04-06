%{
    ===================================================================================
                     Universidad Politécnica Metropolitana de Hidalgo
                             Master in Aerospace Engineering
            
            Filename:   LinearSystems.m
                Date:   Jan 18th, 2022 (Update Feb 10th, 2022)
            
            CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
    ===================================================================================
%}

classdef LinearSystems < handle
    
    properties (Constant, Access = private)
        procedures = ["Gauss-Jordan";
                      "Matriz Adjunta";
                      "Eliminación Gaussiana";
                      "Thomas Llewellyn";
                      "Jacobi";
                      "Gauss-Seidel";
                      "Sobre-Relajación Sucesiva";
                      "Gradiente Conjugado"];
    end % of Constant properties
    
    properties (Access = public)
        method, mCoef, xVector, bVector, degree, tolerance, error, iterMax, iteration, w
    end % of private properties
    
    methods (Access = public)
        %{
            Constructor. Instantiate an object of the class.
        %}
        function app = LinearSystems(), end % of Constructor()
    end % of public methods
    
    methods (Access = public)
        %{
            Set the selected method by user.
            @params Object  $app        Class object.
                    Integer $opt        Selected method label.
        %}
        function setMethod(app,opt)
            app.method = opt;
        end % of setMethod()
        
        %{
            Set the selected method by user.
            @params Object  $app       Class object.
            @return array   $x      
                    double  $error 
                    integer $iter
        %}
        function [x,error,iter] = getResults(app)
            x = '';
            error = 0;
            iter = 0;
            
            if ~isnan(app.error)
                error = app.error;
            end
            
            if ~isnan(app.iteration)
                iter = app.iteration;
            end
            
            for i = 1:app.degree
                x = strvcat(x,sprintf('x%02d = %10.4f',i,app.xVector(i)));
            end
            x = cellstr(x);
        end % of getResults()
        
        function degree = getDegree(app)
            degree = app.degree;
        end
        
        %{
            Captures the system of equations and set as properties.
            @params Object  $app       Class object.
        %}
        function setParams(app,coef,b,tolerance,iterMax)
            app.mCoef = coef;
            app.bVector = b;
            app.tolerance = tolerance;
            app.iterMax = iterMax;
            app.splitSystem();
            app.solve();
        end % of setParams()
        
        %{
            Splits a string to register the system of equations and
            set as properties of this class instanced.
            @params Object  $this       Class object.
        %}
        function isCorrect = splitSystem(app)
            try
                nVector = length(app.bVector);
                temp = zeros(nVector,1);
                
                for i = 1:nVector
                    temp(i) = str2double(app.bVector(i));
                end
                
                app.bVector = temp;
                
                nEquations = length(app.mCoef);
                temp = zeros(nEquations,nEquations);
                
                for i = 1:nEquations
                    aux = split(app.mCoef(i),' ');
                    for j = 1:nEquations
                        temp(i,j) = str2double(aux(j));
                    end
                end
                
                app.mCoef = temp;
                app.degree = nEquations;
                app.xVector = zeros(nEquations,1);
                isCorrect = true;
            catch
                isCorrect = false;
            end
        end % of splitSystem()
        
        %{
            @params Object  $app       Class object.
        %}
        function solve(app)
            app.xVector = zeros(app.degree,1);
            
            switch app.method
                case 'Gauss-Jordan'
                    app.xVector = app.GaussJordanMeth() * app.bVector;
                case 'Matriz Adjunta'
                    app.xVector = app.AdjunctMatrixMeth() * app.bVector;
                case 'Eliminación Gaussiana'
                    app.xVector = app.GaussianElimintationMeth();
                case 'Thomas Llewellyn'
                    app.xVector = app.ThomasMeth();
                case 'Jacobi'
                    app.xVector = app.JacobiMeth();
                case 'Gauss-Seidel'
                    app.xVector = app.GaussSeidelMeth();
                case 'Sobre-Relajación Sucesiva'
                    app.xVector = app.SORMeth();
                case 'Gradiente Conjugado'
                    app.xVector = app.conjGradsMeth();
            end
        end % of solve()
        
        %{
            Implements the Gauss-Jordan method to solve the system of equations.
            @params Object  $app        Class object.
            @return matrix  $mInv       Inverse matrix.
        %}
        function mInv = GaussJordanMeth(app)
            n = app.degree;
            M = [app.mCoef, eye(n)];
            
            for j = 1:n
                M(j,:) = M(j,:)/M(j,j);
                for i = 1:n
                    if i ~= j
                        M(i,:) = M(i,:) - M(i,j)*M(j,:);
                    end
                end
            end
            
            mInv = M(:, n+1:2*n);
        end % of GaussJordanMeth()
        
        %{
            Implements the Adjunct matrix method to solve the system of equations.
            @params Object  $app       Class object.
            @return matrix  $mInv       Inverse matrix.
        %}
        function mInv = AdjunctMatrixMeth(app)
            determinant = det(app.mCoef);
            
            if determinant == 0
                input(sprintf('\n\t\t%s',...
                    '¡La matriz de coeficientes no tiene inversa!'));
                mInv = [];
                return
            end
            
            adjunct = adjoint(app.mCoef);
            mInv = adjunct/determinant;
            
        end % of AdjunctMatrixMeth()
        
        %{
            Implements the Gaussian Elimination method to solve the system of equations.
            @params Object  $app       Class object.
            @return matrix  $xVect     Vector of results.
        %}
        function xVect = GaussianElimintationMeth(app)
            M = [app.mCoef, app.bVector];
            xVect = app.bVector;
            n = app.degree;
            
            for j = 1:n-1
                for i = j+1:n
                    M(i,:) = M(i,:) + M(j,:)*(-M(i,j)/M(j,j));
                end
            end

            for i = n:-1:1
                xVect(i) = M(i,n+1);
                for j = i+1:n
                    xVect(i) = xVect(i) - xVect(j)*M(i,j);
                end
                xVect(i) = xVect(i)/M(i,i);
            end
        end % of GaussianElimintationMeth()
        
        %{
            Implements the Thomas method to solve the system of equations.
            @params Object  $app       Class object.
            @return matrix  $xVect     Vector of results.
        %}
        function xVect = ThomasMeth(app)
            A = app.mCoef;
            L = eye(app.degree);
            U = zeros(app.degree,app.degree);
            
            for i = 2:app.degree
                A(i,i-1) = A(i,i-1)/A(i-1,i-1);
                A(i,i) = A(i,i) - A(i,i-1)*A(i-1,i);
                L(i,i-1) = A(i,i-1);
                U(i-1:i,i) = A(i-1:i,i);
            end
            U(1,1) = A(1,1);
            dVect = L\app.bVector;
            xVect = U\dVect;
        end % of ThomasMeth()
        
        %{
            Implements the Jacobi method to solve the system of equations.
            @params Object  $app       Class object.
            @return matrix  $xVect      Vector of results.
        %}
        function xVect = JacobiMeth(app)
            A = app.mCoef;
            xVect = ones(app.degree,1);
            yVect = xVect;
            app.iteration = 0;
            app.error = 100;
            
            while app.error >= app.tolerance && app.iteration <= app.iterMax
                app.iteration = app.iteration + 1;
                for i = 1:app.degree
                    j = [i:i-1 i+1:app.degree];
                    yVect(i) = (app.bVector(i)-A(i,j)*xVect(j))/A(i,i);
                end
                app.error = max(abs(xVect-yVect))/max(abs(yVect));
                xVect = yVect;
            end
        end % of JacobiMeth()
        
        %{
            Implements the Gauss-Seidel method to solve the system of equations.
            @params Object  $app       Class object.
            @return matrix  $xVect      Vector of results.
        %}
        function xVect = GaussSeidelMeth(app)
            A = app.mCoef;
            xVect = zeros(app.degree,1);
            app.iteration = 0;
            app.error = 100;
            
            while app.error > app.tolerance && app.iteration <= app.iterMax
                app.iteration = app.iteration + 1;
                yVect = xVect;
                for i = 1:app.degree
                    s = A(i,1:i-1) * yVect(1:i-1) + A(i,i+1:app.degree) * yVect(i+1:app.degree);
                    xVect(i) = (app.bVector(i)-s)/A(i,i);
                end
                app.error = norm(yVect-xVect,inf);
            end
        end % of GaussSeidelMeth()
        
        %{
            Implements the Sucessive Over-Relaxtion method to solve the system of equations.
            @params Object  $app       Class object.
            @return matrix  $xVect      Vector of results.
        %}
        function xVect = SORMeth(app)
            xVect = zeros(app.degree,1);
            app.w = 1.25;
            app.error = 100;
            app.iteration = 0;
            
            while app.error >= app.tolerance && app.iteration < app.iterMax
                app.iteration = app.iteration + 1;
                for i = 1:app.degree
                    j = [i:i-1, i+1:app.degree];
                    su = app.bVector(i) - app.mCoef(i,j) * xVect(j);
                    xi = (1-app.w)*xVect(i) + app.w*su/app.mCoef(i,i);
                    app.error = max(abs(app.xVector(i)-xi),app.error);
                    xVect(i) = xi;
                end
                app.error = app.error/max(abs(xVect));
            end
        end % of SORMeth()
        
        %{
            Implements the Conjugate Gradients method to solve the system of equations.
            @params Object  $app       Class object.
            @return matrix  $xVect      Vector of results.
        %}
        function xVect = conjGradsMeth(app)
            app.iteration = 0;
            xVect = zeros(app.degree,1);
            r = app.bVector - app.mCoef*xVect;
            app.error = r'*r;
            app.tolerance = sqrt(eps)*app.error;
            p = r;
            
            while app.error > app.tolerance && app.iteration <= app.iterMax
                app.iteration = app.iteration + 1;
                v = app.mCoef*p;
                alpha = app.error/(p'*v);
                xVect = xVect + alpha*p;
                r = r - alpha*v;
                beta = app.error;
                app.error = r'*r;
                beta = app.error/beta;
                p = r + beta*p;
            end
        end % of conjGradsMeth()
    end % of protected methods
end % of LinearSystems class
