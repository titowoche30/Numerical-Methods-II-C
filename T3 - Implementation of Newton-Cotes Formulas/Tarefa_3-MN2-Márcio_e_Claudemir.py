import math

def int_closed_degree1(a,b,h):
    return (h/2) * (a+b)

def int_open_degree1(a,b,h):
    return ((3*h)/2) * (a+b)

def int_closed_degree2(a, b, c, h):
    return (h/3) * (a + 4*b + c)  

def int_open_degree2(a,b, c,h):
    return (4/3)*h*(2*a - b + 2*c)
                    
def int_closed_degree3(a, b, c, d,h):
    return (3/8)*h * (a + 3*b + 3*c + 1*d)

def int_open_degree3(a,b, c, d, h):
    return (5/24)*h * (11*a + b + c + 11*d)

def int_closed_degree4(a, b, c, d, e, h):
    return ((2*h)/45) * (7*a + 32*b + 12*c + 32*d + 7*e)

def int_open_degree4(a,b, c, d,e, h):
    return (3/5)*h*( (11/2)*a - 7*b + 13*c - 7*d + (11/2)*e)
    

def integrate(f,degree,closed,a,b):
    '''
    Argumentos:
        f - normal or a lambda function
        degree - substitution polynomial degree 
        closed - boolean value, True = closed e False = open
        a e b - integration limits
    '''
    
    if not 1<=degree<=4: 
        print('degree {} not implemented'.format(degree))
        return 
    
    tolerance=10E-6
    result=0
    difference=1
    n=0
    
    while(difference>tolerance):
        l=0
        aux=result
        result=0
        
        if closed: h = ((b-a)/(2**n)) / degree
        else: h = ((b-a)/(2**n)) / (degree +2)
    
        while(l<2**n):
            xi = a+l*((b-a)/2**n)
            
            if closed:
                if degree==4: 
                    result = result + int_closed_degree4(f(xi),f(xi+h),f(xi+2*h),f(xi+3*h),f(xi+4*h),h)
                elif degree==3:
                    result = result + int_closed_degree3(f(xi),f(xi+h),f(xi+2*h),f(xi+3*h),h)
                elif degree==2:
                    result = result + int_closed_degree2(f(xi),f(xi+h),f(xi+2*h),h)
                else:
                    result = result + int_closed_degree1(f(xi),f(xi+h),h)   
            else:
                if degree==4:
                    result = result + int_open_degree4(f(xi+h),f(xi+2*h),f(xi+3*h),f(xi+4*h),f(xi+5*h),h)
                elif degree==3:
                    result = result + int_open_degree3(f(xi+h),f(xi+2*h),f(xi+3*h),f(xi+4*h),h)
                elif degree==2:
                    result = result + int_open_degree2(f(xi+h),f(xi+2*h),f(xi+3*h),h)
                else:
                    result = result + int_open_degree1(f(xi+h),f(xi+2*h),h)  
                    
            l=l+1
            
        difference=abs(result-aux)
        n=n+1
        print('iteration {} = {}'.format(n,result))
        
    return result
    
if __name__ == '__main__':
    #f = lambda x: x**2
    #f = lambda x: 3*x + 7     
    f = lambda x: (math.sin(2*x) + 4*x**2 + 3*x)**2  
    
    degree = int(input('Enter the substitution polynomial degree(1,2,3 or 4)\n'))
    abordagem = str(input('Enter the approach(closed or open)\n'))
    abordagem = abordagem.lower()
    
    if abordagem == 'closed': bol = True
    else: bol = False
  
    result=integrate(f,degree,bol,0,10)
    print('\nresult =',result)