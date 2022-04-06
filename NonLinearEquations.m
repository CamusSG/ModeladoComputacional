%{
    ===================================================================================
                     Universidad Politécnica Metropolitana de Hidalgo
                             Master in Aerospace Engineering
            
            Filename:   NonLinearEquations.m
                Date:   Mar 8th, 2022 (Update Mar 15th, 2022)
            
            CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
    ===================================================================================
%}

classdef NonLinearEquations < handle
    
    properties (Constant, Access = private)
        procedures = ["Bisección";
                      "Falsa Posición";
                      "Punto Fijo";
                      "Newton-Raphson";
                      "Newton-Raphson modificado";
                      "Secante";
                      "Bairstow"];
    end % of Constant properties
    
    properties (Access = public)
        method, equation, interval, tolerance, error, iterMax, iteration, root
    end % of private properties
    
    methods (Access = public)
        %{
            Constructor. Instantiate an object of the class.
        %}
        function app = NonLinearEquations(), end % of Constructor()
    end % of public methods
    
    methods (Access = public)
        %{
            Set the selected method by user.
            @params Object  $app        Class object.
                    Integer $opt        Selected method label.
        %}
        function setMethod(app,tag)
            app.method = tag;
        end % of setMethod()
        
        %{
            Set the selected method by user.
            @params Object  $app       Class object.
            @return double  $root     
                    double  $error 
                    integer $iter
        %}
        function [root,error,iter] = getResults(app)
            root = str2double(sprintf('%.6f',app.root));
            error = str2double(sprintf('%.9f',app.error));
            iter = str2double(sprintf('%d',app.iteration));
        end % of getResults()
        
        %{
            Captures the system of equations and set as properties.
            @params Object  $app       Class object.
        %}
        function setParams(app,equation,interval,tolerance,iterMax)
            app.equation = equation;
            app.interval = interval;
            app.tolerance = tolerance;
            app.iterMax = iterMax;
        end % of setParams()
        
        %{
            @params Object  $app       Class object.
        %}
        function solve(app)
            switch app.method
                case 'Bisección'
                    app.BisectionFalsePosition(1)
                case 'Falsa Posición'
                    app.BisectionFalsePosition(2)
                case 'Punto Fijo'
                    app.FixedPoint()
                case 'Newton-Raphson'
                    app.NewtonRaphson()
                case 'Newton-Raphson Modificado'
                    app.ModifNewtonRaphson()
                case 'Secante'
                    app.Secant()
            end
        end % of solve()
        
        function rst = eval(app,point)
            rst = subs(app.equation,point);
        end
        
        function BisectionFalsePosition(app,opt)
            middle = @(a,b) (a + b)/2;
            calculate = @(a,b) (a*app.eval(b)-(b*app.eval(a)))/(app.eval(b)-app.eval(a));
            lower = app.interval(1);
            upper = app.interval(2);
            app.iteration = 0;
            previous = 0;
            
            while true
                if opt == 1
                    current = middle(lower,upper);
                else
                    current = double(calculate(lower,upper));
                end
                
                if current ~= 0
                    app.error = NonLinearEquations.getError(current, previous);
                else
                    app.error = 0;
                end
                
                app.root = current;

                if app.error < app.tolerance
                    break
                end
                
                if (app.eval(lower)*app.eval(current)) < 0
                    upper = current;
                else
                    lower = current;
                end
                previous = current;
                app.iteration = app.iteration + 1;
            end
        end
        
        function FixedPoint(app)
            app.iteration = 0;
            app.error = 100;
            lower = app.interval(1);
            iter_max = app.iterMax;
            
            while true
                current = double(app.eval(lower));
                app.error = double(NonLinearEquations.getError(current, lower));
                lower = current;
  
                if app.error < app.tolerance || app.iteration >= iter_max
                    app.root = current;
                    break
                end
                app.iteration = app.iteration + 1;
            end
        end
        
        function NewtonRaphson(app)           
            app.iteration = 1;
            app.error = 100;
            func = str2sym('x') - app.equation/diff(app.equation,'x');
            lower = app.interval(1);
            iter_max = app.interval(2);
            
            while true
                current = eval(subs(func,lower));   
                app.error = abs(current - lower);
                lower = current;
                
                if app.error < app.tolerance || app.iteration > iter_max
                    app.root = current;
                    break
                end
                app.iteration = app.iteration + 1;
            end
        end
        
        function ModifNewtonRaphson(app)
            app.iteration = 0;
            app.error = 100;
            f = app.equation;
            df = diff(f,'x');
            ddf = diff(df,'x');
            lower = app.interval(1);
            iter_max = app.interval(2);
            
            func = str2sym('x') - f*df/((df)^2 - f*ddf);
            
            while app.error > app.tolerance || app.iteration < iter_max
                app.iteration = app.iteration + 1;
                current = eval(subs(func,lower));
                app.error = abs(current - lower);
                lower = current;
                app.root = current;
            end
        end
        
        function Secant(app)
            get = @(val) subs(app.equation,val);
            
            app.iteration = 1;
            app.error = 100;
            prev = app.interval(1);
            curr = app.interval(2);
            
            while true
                new = curr-(get(curr)*(prev-curr))/(get(prev)-get(curr));
                app.error = double(abs(new - curr));
                prev = curr;
                curr = new;
                
                if app.error < app.tolerance || app.iteration >= 50
                    app.root = double(curr);
                    break
                end
                app.iteration = app.iteration + 1;
            end
        end
    end % of protected methods
    
    methods (Static)        
        function rst = getError(current,previous)
            rst = abs((current - previous)/current);
        end % of getError()
    end % of Static methods
end % of NonLinearEquations class
