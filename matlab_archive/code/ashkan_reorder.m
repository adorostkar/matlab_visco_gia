function [nlist, flag] = ashkan_reorder(list,Node, np)
    k = 0;
    flag = 0;
    nlist = list;
    tol = 1e-3;
    
    if np == 4
        n1 = Node(:,list(1));
        n2 = Node(:,list(2));
        n3 = Node(:,list(3));
        n4 = Node(:,list(4));
        major = list(:);
        
        
        while(abs(n1(1) - n2(1)) >= tol ||...
                abs(n4(1) - n3(1)) >= tol ||...
                abs(n1(2) - n4(2)) >= tol ||...
                abs(n2(2) - n3(2)) >= tol ||...
                n1(2) > n2(2)             ||...
                n1(1) > n4(1)             &&...
                k < 5)
            
            major = circshift(major,1);
            n1 = Node(:,major(1));
            n2 = Node(:,major(2));
            n3 = Node(:,major(3));
            n4 = Node(:,major(4));
            k = k+1;
        end
        if (k==4)
            flag = 1;
            return;
        end
        nlist = major;
    elseif np == 9
        
        n1 = Node(:,list(1));
        n2 = Node(:,list(2));
        n3 = Node(:,list(3));
        n4 = Node(:,list(4));
        n5 = Node(:,list(5));
        n6 = Node(:,list(6));
        n7 = Node(:,list(7));
        n8 = Node(:,list(8));
        n9 = Node(:,list(9));
        
        major = list(1:4);
        minor = list(5:8);
        
        while(abs(n1(1) - n5(1)) >= tol      &&...
                abs(n1(2) - n5(2)) >= tol && k < 5)
            minor = circshift(minor,1);
            n5 = Node(:,minor(1));
            k = k+1;
        end
        if (k==4)
            flag = 1;
            return;
        end
        k = 0;
        %         [Node(:,major(:)) Node(:,minor(:))]
        while(abs(n1(1) - n5(1)) >= tol ||...
                abs(n1(1) - n2(1)) >= tol ||...
                abs(n8(1) - n9(1)) >= tol ||...
                abs(n8(1) - n6(1)) >= tol ||...
                abs(n4(1) - n7(1)) >= tol ||...
                abs(n4(1) - n3(1)) >= tol ||...
                abs(n1(2) - n8(2)) >= tol ||...
                abs(n1(2) - n4(2)) >= tol ||...
                abs(n5(2) - n9(2)) >= tol ||...
                abs(n5(2) - n7(2)) >= tol ||...
                abs(n2(2) - n6(2)) >= tol ||...
                abs(n2(2) - n3(2)) >= tol ||...
                n1(2) > n5(2)             ||...
                n5(2) > n2(2)             ||...
                n1(2) > n2(2)             ||...
                n1(1) > n8(1)             ||...
                n8(1) > n4(1)             ||...
                n1(1) > n4(1) &&...
                k < 5)
            
            major = circshift(major,1);
            minor = circshift(minor,1);
            n1 = Node(:,major(1));
            n2 = Node(:,major(2));
            n3 = Node(:,major(3));
            n4 = Node(:,major(4));
            n5 = Node(:,minor(1));
            n6 = Node(:,minor(2));
            n7 = Node(:,minor(3));
            n8 = Node(:,minor(4));
            k = k+1;
        end
        if (k==4)
            flag = 1;
            return;
        end
        nlist = [major; minor; list(9)];
    end
end