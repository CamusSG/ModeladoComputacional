%{
    ===================================================================================
                     Universidad Politécnica Metropolitana de Hidalgo
                             Master in Aerospace Engineering
            
            Filename:   DiferentialEquations.m
                Date:   Apr 1st, 2022
            
            CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
    ===================================================================================
%}

classdef DifferentialEquations < handle
    
    properties (Constant, Access = private)
        procedures = ["Euler";
                      "Euler Mejorado";
                      "Runge-Kutta 4to Orden"];
    end % of Constant properties
    
    properties (Access = public)
        method, equation, interval, steps, error, initial, exact, result, ratio
    end % of private properties
    
    methods (Access = public)
        %{
            Constructor. Instantiate an object of the class.
        %}
        function app = DifferentialEquations(), end % of Constructor()
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
            @return double  $result     
                    double  $error 
        %}
        function [result,error,ratio] = getResults(app)
            result = str2double(sprintf('%.16f',app.result));
            error = str2double(sprintf('%.10f',app.error));
            ratio = str2double(sprintf('%.6f',app.ratio));
        end % of getResults()
        
        %{
            Captures the system of equations and set as properties.
            @params Object  $app       Class object.
        %}
        function setParams(app,equation,interval,initial,exact,steps)
            app.equation = equation;
            app.interval = interval;
            app.initial = initial;
            app.exact = exact;
            app.steps = steps;
        end % of setParams()
        
        %{
            @params Object  $app       Class object.
        %}
        function solve(app)
            switch app.method
                case 'Euler'
                    app.Euler()
                case 'Euler Mejorado'
                    app.ImprovedEuler()
                case 'Runge-Kutta 4to Orden'
                    app.RungeKutta4()
            end
        end % of solve()
        
        function Euler(app)
            xPrev = app.interval(1);
            yPrev = app.initial;
            h = (app.interval(2) - app.interval(1)) / app.steps;
            
            for i = 1:app.steps
                f = app.evaluate(xPrev,yPrev);
                
                xCurrent = xPrev + h;
                yCurrent = yPrev + h*f;
                app.error = abs((yCurrent - app.exact)/app.exact) * 100;
                
                xPrev = xCurrent;
                yPrev = yCurrent;
                
                if i == app.steps
                    app.result = yCurrent;
                    app.ratio = h;
                end
            end
        end % of Euler()
        
        function ImprovedEuler(app)
            xPrev = app.interval(1);
            yPrev = app.initial;
            h = (app.interval(2) - app.interval(1)) / app.steps;
            
            for i = 1:app.steps
                f = app.evaluate(xPrev,yPrev);
                xCurrent = xPrev + h;
                zCurrent = yPrev + h*f;
                fCurrent = app.evaluate(xCurrent,zCurrent);                
                yCurrent = yPrev + (h/2)*(f + fCurrent);
               
                app.error = abs((yCurrent - app.exact)/app.exact) * 100;
                
                xPrev = xCurrent;
                yPrev = yCurrent;
                
                if i == app.steps
                    app.result = yCurrent;
                    app.ratio = h;
                end
            end
        end % of ImprovedEuler()
        
        function RungeKutta4(app)
            xPrev = app.interval(1);
            yPrev = app.initial;
            h = (app.interval(2) - app.interval(1)) / app.steps;
            
            for i = 1:app.steps
                p1 = app.evaluate(xPrev,yPrev);
                q1 = app.evaluate(xPrev + 0.5*h,yPrev + 0.5*h*p1);
                r1 = app.evaluate(xPrev + 0.5*h,yPrev + 0.5*h*q1);
                s1 = app.evaluate(xPrev + h,yPrev + h*r1);
                
                xCurrent = xPrev + h;
                yCurrent = yPrev + (h/6)*(p1 + 2*q1 + 2*r1 + s1);
               
                app.error = abs((yCurrent - app.exact)/app.exact) * 100;
                
                xPrev = xCurrent;
                yPrev = yCurrent;
                
                if i == app.steps
                    app.result = yCurrent;
                    app.ratio = h;
                end
            end
        end % of RungeKutta()
        
        function rst = evaluate(app,x,y)
            rst = subs(app.equation,'x',x);
            rst = subs(rst,'y',y);
            rst = eval(rst);
        end
    end
end
