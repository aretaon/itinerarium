# -*- coding: utf-8 -*-
# Generated by Django 1.9.2 on 2016-02-25 02:26
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('itinerarium', '0015_auto_20160225_0324'),
    ]

    operations = [
        migrations.AlterField(
            model_name='post',
            name='image',
            field=models.ImageField(null=True, upload_to='static/images/%Y/%m'),
        ),
    ]
