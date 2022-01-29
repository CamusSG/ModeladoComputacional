%{
    ===================================================================================
                        Universidad Politécnica Metropolitana de Hidalgo
                               Master in Aerospace Engineering

            Filename:   solveLinealSystems.m
                Date:   Jan 18th, 2022 (Update on Jan 28th, 2022)
              Author:   Oscar Federico García Castro
             Contact:   213220003@upmh.edu.mx

            CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
    ===================================================================================
%}

classdef solveLinealSystems < handle
    
    properties (Constant, Access = private)
        procedures = ["Método de Gauss-Jordan";                 "Método de matriz adjunta";
                      "Método de eliminación Gaussiana";        "Método de Thomas Llewellyn";
                      "Método de Jacobi";                       "Método de Gauss-Seidel";
                      "Método de Sobre-Relajación Sucesiva";    "Método del Gradiente Conjugado";
                      "Regresar al menú anterior"];
        actions = ["Usar el mismo sistema de ecuaciones";
                   "Recapturar el sistema de ecuaciones";
                   "Salir del programa"];
    end % of Constant properties
    
    properties (Access = private)
        selMethod, isNew, tagMethod, mCoef, xVector, bVector, degree, mSaved, bSaved, tolerance, error, iterMax, iter, w
    end % of private properties
    
    methods (Access = public)
        %{
            Constructor. Instantiate an object of the class.
        %}
        function this = directMeths()
            this.isNew = true;
            this.init();
        end % of directMeths()
    end % of public methods
    
    methods (Access = protected)
        %{
            Calls functions to capture and solve the system of equations.

            @params Object  $this       Class object.
        %}
        function init(this)
            while true
                if this.isNew
                    this.setEquations();
                else
                    action = directMeths.printMenu('¿QUÉ DESEA HACER?',this.actions);
                    
                    if action == 0
                        break
                    elseif action == 1
                        this.mCoef = this.mSaved;
                        this.bVector = this.bSaved;
                    elseif action == 2
                    	this.setEquations();
                    end
                end
                                
                this.selMethod = directMeths.printMenu('MÉTODOS DE RESOLUCIÓN',this.procedures);
                
                if this.selMethod == 0
                    continue
                elseif this.selMethod >= 5 && this.selMethod <= 8
                    
                    this.iterMax = round(directMeths.getData('Cantidad máxima de interaciones'));
                    
                    if this.selMethod ~= 8
                        this.tolerance = directMeths.getData('Tolerancia');
                    end
                    
                    if this.selMethod == 7
                        this.w = directMeths.getData('Indique el factor de relajación');
                    end
                end
                
                this.solve();
                this.tagMethod = upper(string(this.procedures(this.selMethod)));
                
                if this.selMethod >= 5
                    if this.error > this.tolerance || this.iter == this.iterMax
                        input(sprintf('\n\n%85s.','¡ERROR! No se pudo resolver el sistema de ecuaciones'));
                        continue
                    end
                end
                
                if this.selMethod < 5
                    directMeths.cTitle(sprintf('RESUELTO POR EL %s',this.tagMethod));
                else
                    directMeths.cTitle(sprintf('RESUELTO POR EL %s (ITER. No. %03i, ERROR %6.5f)',this.tagMethod,this.iter,this.error));
                end
                
                directMeths.printSolution(this.degree,this.xVector);
                    
                again = input(sprintf('\n\n%85s: ','¿Desea resolver otro sistema? (Y/N)'),'s');
                
                if strcmp(again,'Y') || strcmp(again,'y')
                    continue
                else
                    break
                end
            end
        end
        
        %{
            Captures the system of equations and set as properties.

            @params Object  $this       Class object.
        %}
        function setEquations(this)
            while true
                directMeths.cTitle('CAPTURANDO EL SISTEMA DE ECUACIONES');
                fprintf('\n\n\t%s\n\n\t\t\t%s\n\n\t%s\n\n\t\t\t%s\n\n', 'Escriba las ecuaciones de la siguiente manera', ...
                    'a11,a12,a13,b1;a21,a22,a23,b2;a31,a32,a33,b3','Por ejemplo','2,3,1,6;3,-2,-4,9;5,-1,-1,4');
                
                try
                    str = input(sprintf('\tCapture el sistema: '),'s');
                catch
                end
                
                if ~this.splitSystem(str)
                    input(sprintf('\n\t\t%s','¡Error al capturar la información! Reintente.'));
                    continue
                end
                
                directMeths.cTitle('SISTEMA DE ECUACIONES CAPTURADA'); 
                directMeths.printSystem(this.mCoef,this.bVector,this.degree);
                
                isCorrect = input(sprintf('\n\n%85s: ','¿La información es correcta? (Y/N)'),'s');
                
                if strcmp(isCorrect,'Y') || strcmp(isCorrect,'y')
                    this.isNew = false;
                    break
                else
                    continue
                end
            end
        end % of setEquations()
        
        %{
            Splits a string to register the system of equations and
            set as properties of this class instanced.

            @params Object  $this       Class object.
        %}
        function isCorrect = splitSystem(this,str)
            try
                str_split = split(str,';');
                items = size(str_split);

                this.degree = items(1);
                this.mCoef = zeros(items(1),items(1));
                this.xVector = zeros(items(1),1);
                this.bVector = zeros(items(1),1);

                for j = 1:items(1)
                    substr = split(str_split(j),',');
                    elements = items(1) + 1;
                    aux = zeros(1,elements);
                    for k = 1:elements(1)
                        aux(k) = str2double(cell2mat(substr(k)));
                    end
                    this.mCoef(j,:) = aux(1:items(1));
                    this.bVector(j) = aux(elements);
                end
                isCorrect = true;
                this.mSaved = this.mCoef;
                this.bSaved = this.bVector;
            catch
                isCorrect = false;
            end
        end % of splitSystem()
        
        %{
            @params Object  $this       Class object.
        %}
        function solve(this)
            switch this.selMethod
                case 1 % Gauss-Jordan Method
                    this.xVector = this.GaussJordanMeth() * this.bVector;
                case 2 % Adjunct Matrix Method
                    this.xVector = this.AdjunctMatrixMeth() * this.bVector;
                case 3 % Gaussian Elimination Method
                    this.xVector = this.GaussianElimintationMeth();
                case 4 % Thomas Llewellyn Method
                    this.xVector = this.ThomasMeth();
                case 5 % Jacobi
                    this.xVector = this.JacobiMeth();
                case 6 % Gauss-Seidel
                    this.xVector = this.GaussSeidelMeth();
                case 7 % Succesive Over-Relaxtion
                    this.xVector = this.SORMeth();
                case 8 % Conjugate Gradients
                    this.xVector = this.conjGradsMeth();
            end
        end % of solve()
        
        %{
            Implements the Gauss-Jordan method to solve the system of equations.

            @params Object  $this       Class object.
            @return matrix  $mInv       Inverse matrix.
        %}
        function mInv = GaussJordanMeth(this)
            n = this.degree;
            M = [this.mCoef, eye(n)];
            
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

            @params Object  $this       Class object.
            @return matrix  $mInv       Inverse matrix.
        %}
        function mInv = AdjunctMatrixMeth(this)
            determinant = det(this.mCoef);
            
            if determinant == 0
                input(sprintf('\n\t\t%s',...
                    '¡La matriz de coeficientes no tiene inversa!'));
                mInv = [];
                return
            end
            
            adjunct = adjoint(this.mCoef);
            mInv = adjunct/determinant;
            
        end % of AdjunctMatrixMeth()
        
        %{
            Implements the Gaussian Elimination method to solve the system of equations.

            @params Object  $this       Class object.
            @return matrix  $xVect      Vector of results.
        %}
        function xVect = GaussianElimintationMeth(this)
            M = [this.mCoef, this.bVector];
            xVect = this.bVector;
            n = this.degree;
            
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

            @params Object  $this       Class object.
            @return matrix  $xVect      Vector of results.
        %}
        function xVect = ThomasMeth(this)
            A = this.mCoef;
            L = eye(this.degree);
            U = zeros(this.degree,this.degree);
            
            for i = 2:this.degree
                A(i,i-1) = A(i,i-1)/A(i-1,i-1);
                A(i,i) = A(i,i) - A(i,i-1)*A(i-1,i);
                L(i,i-1) = A(i,i-1);
                U(i-1:i,i) = A(i-1:i,i);
            end
            U(1,1) = A(1,1);
            dVect = L\this.bVector;
            xVect = U\dVect;
        end % of ThomasMeth()
        
        %{
            Implements the Jacobi method to solve the system of equations.

            @params Object  $this       Class object.
            @return matrix  $xVect      Vector of results.
        %}
        function xVect = JacobiMeth(this)
            A = this.mCoef;
            xVect = ones(this.degree,1);
            this.iter = 0;
            this.error = 100;
            
            while this.error > this.tolerance && this.iter <= this.iterMax
                this.iter = this.iter + 1;
                yVect = xVect;
                for i = 1:this.degree
                    s = A(i,1:i-1) * yVect(1:i-1) + A(i,i+1:this.degree) * yVect(i+1:this.degree);
                    xVect(i) = (this.bVector(i)-s)/A(i,i);
                end
                
                this.error = norm(yVect-xVect,inf);
            end
        end % of JacobiMeth()
        
        %{
            Implements the Gauss-Seidel method to solve the system of equations.

            @params Object  $this       Class object.
            @return matrix  $xVect      Vector of results.
        %}
        function xVect = GaussSeidelMeth(this)
            A = this.mCoef;
            xVect = ones(this.degree,1);
            this.iter = 0;
            this.error = 100;
            
            while this.error > this.tolerance && this.iter <= this.iterMax
                this.iter = this.iter + 1;
                yVect = xVect;
                for i = 1:this.degree
                    s = A(i,1:i-1) * yVect(1:i-1) + A(i,i+1:this.degree) * yVect(i+1:this.degree);
                    xVect(i) = (this.bVector(i)-s)/A(i,i);
                end
                this.error = norm(yVect-xVect,inf);
            end
        end % of GaussSeidelMeth()
        
        %{
            Implements the Sucessive Over-Relaxtion method to solve the system of equations.

            @params Object  $this       Class object.
            @return matrix  $xVect      Vector of results.
        %}
        function xVect = SORMeth(this)
            xVect = zeros(this.degree,1);
            this.error = 100;
            this.iter = 0;
            
            while this.error > this.tolerance && this.iter <= this.iterMax
                this.iter = this.iter + 1;
                for i = 1:this.degree
                    j = [i:i-1, i+1:this.degree];
                    su = this.bVector(i) - this.mCoef(i,j) * xVect(j);
                    xi = (1-this.w)*xVect(i) + this.w*su/this.mCoef(i,i);
                    this.error = max(abs(this.xVector(i)-xi),this.error);
                    xVect(i) = xi;
                end
                this.error = this.error/max(abs(xVect));
            end
        end % of SORMeth()
        
        %{
            Implements the Conjugate Gradients method to solve the system of equations.

            @params Object  $this       Class object.
            @return matrix  $xVect      Vector of results.
        %}
        function xVect = conjGradsMeth(this)
            this.iter = 0;
            xVect = zeros(this.degree,1);
            r = this.bVector - this.mCoef*xVect;
            this.error = r'*r;
            this.tolerance = sqrt(eps)*this.error;
            p = r;
            
            while this.error > this.tolerance && this.iter <= this.iterMax
                this.iter = this.iter + 1;
                v = this.mCoef*p;
                alpha = this.error/(p'*v);
                xVect = xVect + alpha*p;
                r = r - alpha*v;
                beta = this.error;
                this.error = r'*r;
                beta = this.error/beta;
                p = r + beta*p;
            end
        end % of conjGradsMeth()
    end % of protected methods
    
    methods (Static)
        %{
            Displays a menu and requests a choice from the user.

            @params string  $title          Message to be displayed.
                    array   $items          Items that make up the menu.
            @return integer $option         Selected option.
        %}
        function option = printMenu(title,items)
            option = 0;
            nItems = length(items);
    
            while true
                directMeths.cTitle(title);
        
                fprintf('\n\n')
                for i = 1:nItems
                    if i == nItems
                        number = 0;
                    else
                        number = i;
                    end
                    fprintf('\t\t%2i) %s\n',number,items(i))
                end
                directMeths.printBar(90)
                fprintf('\n\n')
    
                try
                    option = input(sprintf('%85s: ','Seleccione una opción'));
            
                    if option >= 0 && option <= nItems-1
                        return
                    else
                        continue
                    end
                catch
                end
            end
        end % of printMenu()
        
        %{
            Print a console-centric title.

            @params string  $title      Message to be displayed.
        %}
        function cTitle(title)
            clc;
            str = strjust(sprintf('%86s',title),'center');
            directMeths.printBar(90)
            fprintf('\n\t||%s||',str)
            directMeths.printBar(90)
        end % of cTitle()
        
        %{  
            Prints a bar with equal symbols.

            @params integer  $line      Number of items in the bar.
        %}
        function printBar(line)
            fprintf('\n\t')
            for i = 1:line
                fprintf('=')
            end
        end % of printBar() 
        
        %{
            Prints a system of equations of degree 'n'.

            @params matrix  $A      Matrix of coefficients
                    array   $b      Vector of results
                    integer $n      Degree of the system
        %}
        function printSystem(A,b,n)
            disp(' ')
            for i = 1:n
                for j = 1:n
                    if j == 1
                        str = sprintf('\n\t\t|%10.3f ',A(i,j));
                    elseif j == n
                        str = sprintf('%10.3f|\t|x%02i|\t=\t|%10.3f|',A(i,j),i,b(i));
                    else
                        str = sprintf('%10.3f ',A(i,j));
                    end
                    fprintf('%s',str);
                end
            end
            disp(' ')
            directMeths.printBar(90);
        end % of printSystem()
        
        %{
            @params integer $n      Degree of the system
                    array   $x      Vector of results
        %}
        function printSolution(n,x)
            disp(' ')
            for i = 1:n
                str = sprintf('\n%36s|x%02i|\t=\t|%12.5f|',' ',i,x(i));
                fprintf('%s',str);
            end
            disp(' ')
            directMeths.printBar(90);
        end % of showSolution()
        
        %{
            @params string  $msg      Degree of the system
                    data    $double   Vector of results
        %}
        function data = getData(msg)
            while true
                try
                    data= input(sprintf('%85s: ',msg));
            
                    if data > 0
                        return
                    else
                        input(sprintf('\n\t\t%s','¡Error al capturar la información! Reintente.'));
                        continue
                    end
                catch
                end
            end
        end % of getData()
    end % of Static methods
end % of classdef
