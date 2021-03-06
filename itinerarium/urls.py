from django.conf.urls import url
from . import views
from django.contrib.auth import views as auth_views

urlpatterns = [
    url(r'^$', views.post_list, name='post_list'),
    url(r'post/(?P<pk>[0-9]+)/$', views.post_detail, name='post_detail'),	#look for url starting with post/1232/ and transfer 1232 to the variable pk that is then returned to views.post_detail
    url(r'^post/new/$', views.post_new, name='post_new'),
    url(r'^post/(?P<pk>[0-9]+)/edit/$', views.post_edit, name='post_edit'),
    url(r'^drafts/$', views.post_draft_list, name='post_draft_list'),
    url(r'^post/(?P<pk>[0-9]+)/publish/$', views.post_publish, name='post_publish'),
    url(r'^post/(?P<pk>[0-9]+)/remove/$', views.post_remove, name="post_remove"),
    url(r'^post/(?P<pk>[0-9]+)/comment/$', views.add_comment_to_post, name='add_comment_to_post'),
    url(r'^accounts/login/$', auth_views.login),
    url(r'^accounts/logout/$', auth_views.logout, {'next_page': '/'}),
    url(r'^comment/(?P<pk>[0-9]+)/approve/$', views.comment_approve, name="comment_approve"),
    url(r'^comment/(?P<pk>[0-9]+)/remove/$', views.comment_remove, name="comment_remove"),
]
